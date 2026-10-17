#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated Electron Microscopy (EM) Particle Analysis Toolkit

This module provides the EM_Particle_Finder class for particle detection, selection of X-ray spectra (EDS, WDS) acquisition spots,
and particle size/statistics analysis.
It is designed to interface with an EM driver and an EM_Controller object to streamline the workflow from image acquisition to data export.

Currently only supports SEM.

Main Class
----------
EM_Particle_Finder
    Automates particle detection, mask generation, and X-ray spectra spot selection on detected particles.
    Supports both fully automated and manual (user-guided) collection modes.
    Provides methods for collecting particle statistics and managing frame navigation.

Example Usage
-------------
# EDS spot acquisition workflow:
>>> EM_controller = EM_Controller(...)
>>> EM_controller.initialise_SEM()
>>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, ...)
>>> while particle_finder.go_to_next_particle():
...     xy_spot_list = particle_finder.get_XS_acquisition_spots_coord_list()
...     for (x,y) in xy_spot_list:
...         EM_controller.acquire_XS_spectrum()


# Particle analysis workflow:
>>> EM_controller = EM_Controller(...)
>>> EM_controller.initialise_SEM()
>>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, results_dir="./results")
>>> particle_finder.get_particle_stats(n_par_target=500)

Notes
-----
- Requires an initialised EM_Controller object, unless used in development mode
- Requires a working microscope driver and appropriate hardware configuration (see the documentation page
  "Add new Electron Microscope Driver"), in particular: is_at_EM (True if running at the electron microscope),
  get_image_data(width, height, channel) (acquire an image), get_frame_width() (current frame width, in mm) and
  move_to(x, y) (move the stage).

Created on Wed Jul 31 09:28:07 2024

@author: Andrea
"""
# Standard library imports
import os
import time
import random
import warnings

# Third-party imports
import cv2
import numpy as np

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from autoemx.core.em_runtime.image_utilities import particle_mark_annotations
from autoemx.core.em_runtime.xsp_spot_selection import (
    XSpSpotSelectionContext,
    XSpSpotSelectorCallback,
    validate_xsp_spot_pixels,
)

# Local project imports
from autoemx.utils.helper import (
    Prompt_User,
    draw_scalebar,
    print_double_separator,
)
from autoemx.config import (
    PowderMeasurementConfig,
)
from autoemx.config.ledger_schemas import (
    Coordinate2D,
    FrameRecord,
    FrameSegmentation,
    ParticleInfo,
    ParticleSizeBand,
    ParticleStatsFrame,
    SegmentationConfig,
    SizeBand,
    SizeMeasurement,
)
from autoemx.core import particle_statistics
import autoemx.utils.constants as cnst
from autoemx import microscope_drivers as EM_driver
import autoemx.core.em_runtime.particle_segmentation_models as par_seg_models

from autoemx._logging import get_logger
logger = get_logger(__name__)

#: Particles whose bounding box is within this many pixels of the frame edge are left out of particle size statistics
#: (their size is not known), and the others corrected accordingly (see particle_statistics.edge_correction_factor).
FRAME_EDGE_MARGIN_PX = 3

#: A particle size statistics measurement in progress is saved to the sample ledger every this many frames (and at
#: the end of each size band), so that it can be resumed if interrupted.
STATS_SAVE_EVERY_N_FRAMES = 5


def min_valid_diameter_um(pixel_size_um: float, min_par_pixels: int) -> float:
    """Equivalent diameter (µm) of a particle covering ``min_par_pixels`` pixels of size ``pixel_size_um``."""
    return float(2.0 * pixel_size_um * np.sqrt(min_par_pixels / np.pi))


def plan_size_bands_um(
    d_min_um: float,
    d_max_um: float,
    im_width_px: int,
    min_par_pixels: int,
    frame_width_range_um: Optional[Tuple[float, float]] = None,
) -> List[Tuple[float, float, float]]:
    """
    Size bands of a particle size statistics measurement, as ``(d_min, d_max, frame_width)`` in µm, largest first.

    Bands are planned from the smallest particles up, since they are the hardest to find: the first band is
    searched with the largest frame in which particles of ``d_min_um`` still cover ``min_par_pixels`` pixels, and
    goes up to the largest particles fitting that frame (a tenth of its width); the next band starts there, and so
    on. A band whose frame reaches :data:`MAX_FRAME_WIDTH_UM` (or the largest frame of the microscope) takes all the
    particles left, up to ``d_max_um``. Any narrow remainder is thus measured with the largest frames.

    Parameters
    ----------
    d_min_um, d_max_um : float
        Diameter range, in µm.
    im_width_px : int
        Width of the frame images, in pixels.
    min_par_pixels : int
        Minimum number of pixels of a particle for its area to be measured.
    frame_width_range_um : (float, float), optional
        Smallest and largest frame width of the microscope, in µm. If its smallest frame is too large for the
        smallest particles to cover ``min_par_pixels`` pixels, the first band starts above ``d_min_um``.

    Raises
    ------
    ValueError
        If the limits are not finite and positive, or the frames are too narrow in pixels for particles of
        ``min_par_pixels`` pixels to fit them.
    """
    if not np.isfinite(d_min_um) or not np.isfinite(d_max_um):
        raise ValueError("Diameter limits must be finite")
    if d_min_um <= 0 or d_max_um <= 0:
        raise ValueError("Diameter limits must be positive")
    if d_max_um < d_min_um:
        raise ValueError("d_max_um must be greater than or equal to d_min_um")
    # Smallest valid diameter per µm of frame width
    pixel_limit_per_fw = min_valid_diameter_um(1.0 / im_width_px, min_par_pixels)
    if pixel_limit_per_fw * FRAME_WIDTH_PER_MAX_DIAMETER >= 1:
        raise ValueError(
            f"Frames of {im_width_px} pixels are too narrow to measure particles of {min_par_pixels} pixels: "
            "reduce min_par_pixels"
        )
    min_fw_um, max_fw_um = frame_width_range_um if frame_width_range_um is not None else (0.0, np.inf)
    max_fw_um = min(max_fw_um, MAX_FRAME_WIDTH_UM)

    bands: List[Tuple[float, float, float]] = []
    d_lo = float(d_min_um)
    while True:
        fw_um = d_lo / pixel_limit_per_fw  # largest frame for the smallest particles of the band
        is_largest_frame = fw_um >= max_fw_um
        fw_um = max(min(fw_um, max_fw_um), min_fw_um)
        d_lo = max(d_lo, pixel_limit_per_fw * fw_um)
        d_hi = fw_um / FRAME_WIDTH_PER_MAX_DIAMETER
        if is_largest_frame or d_hi >= d_max_um:
            d_hi = float(d_max_um)
        if d_lo > d_hi * (1 + 1e-12):  # only if the microscope cannot zoom in enough for the particles left
            break
        bands.append((d_lo, d_hi, fw_um))
        if d_hi >= d_max_um:
            break
        d_lo = d_hi
    return bands[::-1]


def diameter_um_to_area_um2(diameter_um: float) -> float:
    """Convert equivalent circular diameter (µm) to area (µm²)."""
    return float(np.pi * (diameter_um / 2.0) ** 2)


#: Frame width of particle size statistics relative to the largest particles it measures, and largest frame width.
FRAME_WIDTH_PER_MAX_DIAMETER = 10.0
MAX_FRAME_WIDTH_UM = 500.0


def frame_width_um_from_max_area(max_area_par: float) -> float:
    """Search-frame width (µm) from maximum accepted particle area (µm²)."""
    max_par_radius = np.sqrt(max_area_par / np.pi)
    return float(min(FRAME_WIDTH_PER_MAX_DIAMETER * 2.0 * max_par_radius, MAX_FRAME_WIDTH_UM))


def is_binary_mask(par_mask: np.ndarray) -> bool:
    """Whether a particle mask is binary (bool, or uint8 with values 0 and 255 only), not a label image."""
    if par_mask.dtype == np.bool_:
        return True
    return par_mask.dtype == np.uint8 and bool(np.isin(np.unique(par_mask), (0, 255)).all())


def label_image_stats(labels: np.ndarray) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    Number of labels (background included), statistics and centroids of a label image, as returned by
    cv2.connectedComponentsWithStats (labels absent from the image have zero statistics).
    """
    from scipy import ndimage

    num_labels = int(labels.max()) + 1 if labels.size else 1
    stats = np.zeros((num_labels, 5), dtype=np.int32)
    centroids = np.zeros((num_labels, 2), dtype=np.float64)
    flat = labels.ravel()
    areas = np.bincount(flat, minlength=num_labels)
    ys, xs = np.indices(labels.shape)
    sum_x = np.bincount(flat, weights=xs.ravel(), minlength=num_labels)
    sum_y = np.bincount(flat, weights=ys.ravel(), minlength=num_labels)
    present = areas > 0
    centroids[present, 0] = sum_x[present] / areas[present]
    centroids[present, 1] = sum_y[present] / areas[present]
    stats[:, cv2.CC_STAT_AREA] = areas
    # Bounding boxes: of the particles (labels 1..), and of the background
    slices = ndimage.find_objects(labels)
    for label, slc in enumerate(slices, start=1):
        if slc is None:
            continue
        stats[label, cv2.CC_STAT_TOP], stats[label, cv2.CC_STAT_LEFT] = slc[0].start, slc[1].start
        stats[label, cv2.CC_STAT_HEIGHT] = slc[0].stop - slc[0].start
        stats[label, cv2.CC_STAT_WIDTH] = slc[1].stop - slc[1].start
    if present[0]:
        bg_ys, bg_xs = np.nonzero(labels == 0)
        stats[0, :4] = (bg_xs.min(), bg_ys.min(), bg_xs.max() - bg_xs.min() + 1, bg_ys.max() - bg_ys.min() + 1)
    return num_labels, stats, centroids

#%% Electron Microscope Particle Finder class    
@dataclass
class StatsProgress:
    """
    Particle size statistics measurement being acquired: target, settings, working records of its size bands (with
    the particles of their frames), and its record in the attached ledger with the segmentation used.
    """

    n_par_target: int
    powder_meas_cfg: PowderMeasurementConfig
    bands: List[ParticleSizeBand]
    measurement: Optional[SizeMeasurement] = None
    segmentation_id: Optional[int] = None


class EM_Particle_Finder:
    """
    Automated particle analysis and X-ray spectra (EDS, WDS) acquisition in an electron microscope (EM).

    This class provides methods for particle detection, selection of X-ray spectra (EDS, WDS) acquisition spots,
    and particle size/statistics analysis.
    It is designed to interface with an EM driver and an EM_Controller object to streamline the workflow from
    image acquisition to data export.


    Main Methods
    ------------
    initialise_SEM()
        Wakes up the SEM microscope, sets measurement parameters, and evaluates locations to scan for particles.
    go_to_next_particle()
        Locates and moves to the next detected particle.
    get_XS_acquisition_spots_coord_list()
        Determines X-ray spectra (EDS, WDS) acquisition spots on the currently detected particle.
    get_particle_stats()
        Scans the full sample, collecting particle size distribution statistics.

    Examples
    --------
    # EDS spot acquisition workflow:
    >>> EM_controller = EM_Controller(...)
    >>> EM_controller.initialise_SEM()
    >>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, ...)
    >>> while particle_finder.go_to_next_particle():
    ...     xy_spot_list = particle_finder.get_XS_acquisition_spots_coord_list()
    ...     for (x,y) in xy_spot_list:
    ...         EM_controller.acquire_XS_spectrum()
    
    
    # Particle analysis workflow:
    >>> EM_controller = EM_Controller(...)
    >>> EM_controller.initialise_SEM()
    >>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, results_dir="./results")
    >>> particle_finder.get_particle_stats(n_par_target=500)

    Attributes
    ----------
    EM : EM_Controller
        Reference to the parent EM_Controller instance (must be initialised before use).
    powder_meas_cfg : PowderMeasurementConfig
        Configuration object for powder measurement (see dataclass for details).
    par_selection_mode : str
        Particle navigation mode from ``powder_meas_cfg`` ('auto' or 'manual').
    results_dir : str or None
        Directory for saving result images and data.
    verbose : bool
        If True, print progress and information to the console.
    development_mode : bool
        If True, enables extra visualizations and debug image saving.

    Internal Attributes
    -------------------
    _sample_id : str
        Identifier for the sample. Inherited from EM_Controller.
    _im_width : int
        Image width in pixels. Inherited from EM_Controller.
    _im_height : int
        Image height in pixels. Inherited from EM_Controller.
    tot_par_cntr : int
        Counter for the total number of particles analyzed.
    analyzed_pars : list of ParticleInfo
        Particle geometry records collected during particle size statistics.
        When a ledger is attached, this list is the same object as ``ledger.particles``.

    Notes
    -----
    - This class assumes that an EM_Controller object is already initialised and passed as `EM_controller`.
    - All configuration validation is performed in the respective dataclasses.
    - Attributes with a leading underscore are for internal use and should not be accessed directly by users.
    """
    def __init__(
        self,
        EM_controller,
        powder_meas_cfg: PowderMeasurementConfig,
        xsp_spot_selector: Optional[XSpSpotSelectorCallback] = None,
        results_dir: Optional[str] = None,
        verbose: bool = True,
        development_mode: bool = True
    ):
        """
        Initialize an EM_Particle_Finder object for automated/manual electron microscopy particle analysis.

        This class requires that an EM_Controller object has already been fully initialised and passed as `EM_controller`.

        Parameters
        ----------
        EM_controller : EM_Controller
            Reference to the parent EM_Controller instance (must be initialised before use).
        powder_meas_cfg : PowderMeasurementConfig
            Configuration object for powder measurement (see dataclass for details).
        xsp_spot_selector : callable, optional
            Required when ``powder_meas_cfg.par_spot_selection_mode='callback'``.
        results_dir : str, optional
            Directory to save result images and data (default: None).
        verbose : bool, optional
            If True, print progress and information to the console (default: True).
        development_mode : bool, optional
            If True, enables extra visualizations and debug image saving (default: True).
            
        Raises
        ------
        None

        Notes
        -----
        - All configuration validation is performed in the respective dataclasses.
        - This initializer assumes that configuration dataclasses are valid and complete.
        - Additional internal attributes are initialized for particle tracking and statistics.
        """
        # Electron microscope controller object
        self.EM = EM_controller
        self.powder_meas_cfg = powder_meas_cfg
        
        ### Inherit attributes
        if EM_controller is not None:
            # EM_controller is set to None when simply processing data
            self._sample_id = EM_controller.sample_id
            self._im_width = EM_controller.im_width
            self._im_height = EM_controller.im_height
        
        self.xsp_spot_selector = xsp_spot_selector
        if self.powder_meas_cfg.par_spot_selection_mode == 'callback' and self.xsp_spot_selector is None:
            raise ValueError(
                "powder_meas_cfg.par_spot_selection_mode='callback' requires xsp_spot_selector."
            )
        # NOTE: self.tot_par_cntr is initialized below, so _select_par_prompt_title will be set on first use
        if self.is_manual_particle_selection:
            self._select_par_prompt_title = f"Select position for particle #{self.tot_par_cntr}"
            self._select_par_prompt_message = (
                "Center the image around the particle you want\n"
                "to analyse, then click OK."
            )

        # --- General options
        self.results_dir = results_dir
        self.verbose = verbose
        self.development_mode = development_mode  # Set whether in code development mode. Triggers visualisations and image saving

        # --- Variable initializations
        self.tot_par_cntr = 0  # Keeps track of total number of particles analysed
        self.ref_image = None
        self.current_par_area_um2: Optional[float] = None
        self.current_par_coordinates_mm: Optional[tuple] = None
        self._fr_par_cntr = 0
        self._num_par_in_frame = 0
        self.analyzed_pars: List[ParticleInfo] = []
        self._ledger = None  # Optional SampleLedger for particle-stats / list selection
        self._ledger_path: Optional[str] = None
        self._list_particle_ids: Optional[List[int]] = None
        self._list_particle_coords_mm: Optional[List[Tuple[float, float]]] = None
        self._list_target_index: int = 0
        self._list_mode_absolute_ids: bool = False
        self._particles_by_id: Dict[int, ParticleInfo] = {}
        # Particle size statistics measurement in progress, and the size band being collected
        self._stats_progress: Optional[StatsProgress] = None
        self._stats_band: Optional[ParticleSizeBand] = None
        self._stats_max_n_frames: Optional[int] = None  # frame limit of a measurement stopping on precision
        # Frames searched by EDS counted by the measurement in progress (get_particle_stats(count_eds_frames=True)),
        # and the diameters of the particles counted in them
        self._eds_frames: List[particle_statistics.EdsSearchFrame] = []
        self._eds_diameters_um: np.ndarray = np.array([])
        # Frames searched for EDS (automatic particle selection) not saved to the ledger yet, with the size range
        # they count: ((d_min_um, d_max_um, frame_width_um), frame record)
        # (diameters searched and frame width, frame, outlines of its particles or None)
        self.pending_eds_frames: List[Tuple[Tuple[float, float, float], ParticleStatsFrame, Optional[dict]]] = []

    @property
    def is_manual_particle_selection(self) -> bool:
        """Whether particle navigation uses manual centering prompts."""
        return self.powder_meas_cfg.par_selection_mode == 'manual'

    @property
    def is_list_particle_selection(self) -> bool:
        """Whether particle navigation follows a caller-supplied id/coordinate list."""
        return self.powder_meas_cfg.par_selection_mode == 'list'

    @property
    def uses_absolute_particle_ids(self) -> bool:
        """True when ``tot_par_cntr`` is a ledger particle id (list-id mode)."""
        return bool(self._list_mode_absolute_ids)

    @property
    def par_selection_mode(self) -> str:
        return self.powder_meas_cfg.par_selection_mode

    @staticmethod
    def _normalize_stage_coordinates_mm(
        coordinates: Sequence[float],
    ) -> Tuple[float, float]:
        if coordinates is None or len(coordinates) != 2:
            raise ValueError("Stage coordinates must be a sequence of two finite floats (x, y) in mm")
        x, y = float(coordinates[0]), float(coordinates[1])
        if not (np.isfinite(x) and np.isfinite(y)):
            raise ValueError("Stage coordinates must be finite")
        return (x, y)

    def set_particle_targets(
        self,
        particle_ids: Optional[Iterable[int]] = None,
        particle_coords_mm: Optional[Iterable[Sequence[float]]] = None,
        ledger=None,
    ) -> None:
        """Configure list-mode targets (ledger particle ids XOR stage coordinates).

        Parameters
        ----------
        particle_ids :
            Existing ``ParticleInfo.id`` values. Requires ``ledger``; each id must
            exist and have absolute stage ``coordinates``.
        particle_coords_mm :
            Absolute stage positions ``(x, y)`` in mm for new particles.
        ledger :
            Sample ledger used to resolve and validate ``particle_ids``.
        """
        has_ids = particle_ids is not None
        has_coords = particle_coords_mm is not None
        if has_ids == has_coords:
            raise ValueError(
                "set_particle_targets requires exactly one of particle_ids or "
                "particle_coords_mm"
            )

        self._list_target_index = 0
        self._list_particle_ids = None
        self._list_particle_coords_mm = None
        self._list_mode_absolute_ids = False
        self._particles_by_id = {}

        if has_ids:
            if ledger is None:
                raise ValueError("particle_ids list mode requires a SampleLedger")
            ids = [int(pid) for pid in particle_ids]
            if not ids:
                raise ValueError("particle_ids must be a non-empty list")
            by_id = {int(p.id): p for p in ledger.particles}
            missing = [pid for pid in ids if pid not in by_id]
            if missing:
                raise ValueError(
                    f"Particle id(s) not found in ledger: {missing}"
                )
            missing_coords = [
                pid for pid in ids if by_id[pid].coordinates is None
            ]
            if missing_coords:
                raise ValueError(
                    "Particle id(s) in ledger are missing absolute stage "
                    f"coordinates: {missing_coords}"
                )
            self._list_particle_ids = ids
            self._list_mode_absolute_ids = True
            self._particles_by_id = by_id
            self._ledger = ledger
        else:
            coords = [
                self._normalize_stage_coordinates_mm(c) for c in particle_coords_mm
            ]
            if not coords:
                raise ValueError("particle_coords_mm must be a non-empty list")
            self._list_particle_coords_mm = coords
            self._list_mode_absolute_ids = False
            if ledger is not None:
                self._ledger = ledger

    def _default_particle_frame_width_mm(self) -> float:
        """Zoom FOV (mm) when particle size is unknown (coordinate list mode)."""
        max_radius_um = float(np.sqrt(self.powder_meas_cfg.max_area_par / np.pi))
        fw_um = max_radius_um * 2.0 * 1.8 + 10.0
        return float(fw_um / 1000.0)

    def _frame_width_mm_for_particle(self, particle: ParticleInfo) -> float:
        """Zoom FOV (mm) from an existing ParticleInfo size when available."""
        if particle.eq_diameter_um is not None and particle.eq_diameter_um > 0:
            eq_d_um = float(particle.eq_diameter_um)
        elif particle.area_um is not None and particle.area_um > 0:
            eq_d_um = float(2.0 * np.sqrt(float(particle.area_um) / np.pi))
        else:
            return self._default_particle_frame_width_mm()
        fw_um = eq_d_um * 1.8 + 10.0
        return float(fw_um / 1000.0)

    def _go_to_listed_particle(
        self,
        particle_id: Optional[int] = None,
        coordinates: Optional[Sequence[float]] = None,
    ) -> bool:
        """Move to a listed particle id or absolute stage coordinate."""
        self._check_EM_controller_initialization()

        if particle_id is not None and coordinates is not None:
            raise ValueError("Pass particle_id or coordinates, not both")

        if particle_id is None and coordinates is None:
            if self._list_particle_ids is not None:
                if self._list_target_index >= len(self._list_particle_ids):
                    return False
                particle_id = self._list_particle_ids[self._list_target_index]
                self._list_target_index += 1
            elif self._list_particle_coords_mm is not None:
                if self._list_target_index >= len(self._list_particle_coords_mm):
                    return False
                coordinates = self._list_particle_coords_mm[self._list_target_index]
                self._list_target_index += 1
            else:
                raise RuntimeError(
                    "list particle selection requires set_particle_targets(...) "
                    "before go_to_next_particle()"
                )

        if particle_id is not None:
            pid = int(particle_id)
            particle = self._particles_by_id.get(pid)
            if particle is None and self._ledger is not None:
                particle = next(
                    (p for p in self._ledger.particles if int(p.id) == pid),
                    None,
                )
            if particle is None:
                raise ValueError(f"Particle id {pid} not found in ledger")
            if particle.coordinates is None:
                raise ValueError(
                    f"Particle id {pid} has no absolute stage coordinates in the ledger"
                )
            coords = (
                float(particle.coordinates.x),
                float(particle.coordinates.y),
            )
            frame_width_mm = self._frame_width_mm_for_particle(particle)
            self.tot_par_cntr = pid
            self._list_mode_absolute_ids = True
            self.current_par_area_um2 = (
                float(particle.area_um) if particle.area_um is not None else None
            )
            self.current_par_coordinates_mm = coords
            # Keep ledger ParticleInfo in sync with the navigation target.
            particle.coordinates = Coordinate2D(x=coords[0], y=coords[1])
        else:
            coords = self._normalize_stage_coordinates_mm(coordinates)
            frame_width_mm = self._default_particle_frame_width_mm()
            self.tot_par_cntr += 1
            self._list_mode_absolute_ids = False
            self.current_par_coordinates_mm = coords

        self.EM.move_to_pos(coords)
        self.EM.set_frame_width(frame_width_mm)
        self.EM.pixel_size_um = frame_width_mm / self._im_width * 10**3
        self.EM.adjust_BCF()
        return True

    
    def _check_EM_controller_initialization(self) -> None:
        """
        Check whether the associated EM_Controller instance is initialized.
    
        Raises
        ------
        RuntimeError
            If the EM_Controller is not initialized.
    
        Notes
        -----
        This method should be called before performing any operation that
        requires an initialized microscope.
        """
        if not self.EM.is_initialized:
            raise RuntimeError(
                "EM_Controller is not initialized. Please call initialise_SEM() "
                "(or initialise_STEM(), if supported) before using this method."
            )
            
            
    #%% Particle Navigation
    # =============================================================================
    def go_to_next_particle(
        self,
        particle_id: Optional[int] = None,
        coordinates: Optional[Sequence[float]] = None,
    ):
        '''
        Moves the microscope to the next particle, centering and zooming on it.
        It also re-adjusts brightness, contrast and focus.

        Modes
        -----
        - ``auto``: find the next particle in the search-frame grid.
        - ``manual``: prompt the user to center a particle.
        - ``list``: consume the next target from ``set_particle_targets`` (ledger
          particle id or absolute stage coordinates in mm). Explicit ``particle_id``
          / ``coordinates`` arguments also select a specific target directly.

        Parameters
        ----------
        particle_id : int, optional
            Existing ledger particle id to revisit. Updates that ParticleInfo.
        coordinates : sequence of float, optional
            Absolute stage ``(x, y)`` in mm. Creates a new particle on acquisition.

        Returns
        -------
        bool
            True if it could successfully move to a particle.
            False if no more particles are present in the sample or if execution is stopped by the user.
        '''
        if particle_id is not None and coordinates is not None:
            raise ValueError("Pass particle_id or coordinates, not both")

        # Explicit target or configured list mode.
        if (
            self.is_list_particle_selection
            or particle_id is not None
            or coordinates is not None
        ):
            self.current_par_area_um2 = None
            self.current_par_coordinates_mm = None
            return self._go_to_listed_particle(
                particle_id=particle_id,
                coordinates=coordinates,
            )

        self.current_par_area_um2 = None
        self.current_par_coordinates_mm = None
        if not self.is_manual_particle_selection:
            self._check_EM_controller_initialization()
            
            # Automated search of particles
            # Check if all (or the max allowed amount) particles have been analysed in the current frames.
            # If yes, re-calculate particle positions in the next frame
            if self._fr_par_cntr >= self._num_par_in_frame or self._fr_par_cntr >= self.powder_meas_cfg.max_n_par_per_frame:
                # Loop until a frame with particles is found
                were_particles_found = False
                while not were_particles_found:
                    were_particles_found = self._get_particles_coordinates_in_frame()
                    # Check if all frames in sample have been analysed
                    if were_particles_found is None:
                        return False
    
            # Move center of image to the centroid of the particle
            par_pos = self.par_pos_abs_mm[self._fr_par_cntr]
            self.EM.move_to_pos(par_pos)
            self.current_par_coordinates_mm = (float(par_pos[0]), float(par_pos[1]))
    
            # Set frame width to zoom around particle
            frame_width_mm = self.par_fw_mm[self._fr_par_cntr]
            self.EM.set_frame_width(frame_width_mm)
    
            # Adjust EM settings (focus, contrast, brightness)
            self.EM.adjust_BCF()
    
            # Update counter of particles in current frame
            self._fr_par_cntr += 1
    
            # Update counter of total number of particles analysed
            self.tot_par_cntr += 1
    
        else:  # Manually look for particle
            prompt = Prompt_User(self._select_par_prompt_title, self._select_par_prompt_message)
            prompt.run()
    
            if prompt.execution_stopped:  # Check if execution was stopped after the loop
                logger.warning("⚠️ Execution stopped by the user.")
                return False
    
            if prompt.ok_pressed:
                frame_width_mm = EM_driver.get_frame_width()
                self.EM.pixel_size_um = frame_width_mm / self._im_width * 10**3  # um
                self.tot_par_cntr += 1
                current_pos = getattr(getattr(self.EM, "frame_navigator", None), "_current_pos", None)
                if current_pos is not None:
                    self.current_par_coordinates_mm = (float(current_pos[0]), float(current_pos[1]))
    
        return True
    
    
    def _stats_frame_id(self) -> Optional[str]:
        """
        Id of the current frame: its label, made unique (``#2``, ``#3``, ...) in the attached ledger (or else in the
        size band being collected) if a frame with the same label was recorded (e.g. by an earlier measurement).
        """
        label = self.EM.current_frame_label
        if label is None:
            return None
        frame_id = str(label)
        if self._ledger is not None:
            return self._ledger.unique_frame_id(frame_id)
        band = self._stats_band
        if band is not None:
            recorded = {frame.frame_id for frame in band.frames}
            n = 2
            unique_id = frame_id
            while unique_id in recorded:
                unique_id, n = f"{frame_id}#{n}", n + 1
            frame_id = unique_id
        return frame_id

    def _unique_image_stem(self, stem: str) -> str:
        """Name of the images of a frame (``<stem>_particles``), made unique among the saved images (``#2``, ...)."""
        folder = self.results_dir
        if not folder or not os.path.isdir(folder):
            return stem
        names = {os.path.splitext(name)[0] for name in os.listdir(folder)}
        unique, n = stem, 2
        while f"{unique}_particles" in names:
            unique, n = f"{stem}#{n}", n + 1
        return unique

    def _valid_particles(self, frame_image: np.ndarray, stats: np.ndarray, centroids: np.ndarray, num_labels: int):
        """
        Particles of a frame valid for particle size statistics (within the area limits, of at least min_par_pixels
        pixels, not touching the frame edge): their areas (µm²), frame-edge corrections and marks (x, y, radius in
        pixels), and their connected-component labels.
        """
        im_height, im_width = frame_image.shape[:2]
        pixel_size_um = float(self.EM.pixel_size_um)
        areas_um2, corrections, marks, labels = [], [], [], []
        for i in range(1, num_labels):  # Skip the background component (index 0)
            area_px = int(stats[i, cv2.CC_STAT_AREA])
            if (
                area_px >= self.powder_meas_cfg.min_par_pixels  # enough pixels for a valid area
                and self._is_particle_area_ok(area_px)
                and not self.is_particle_at_frame_edge(stats, i)
            ):
                areas_um2.append(float(area_px * pixel_size_um**2))
                corrections.append(particle_statistics.edge_correction_factor(
                    int(stats[i, cv2.CC_STAT_WIDTH]), int(stats[i, cv2.CC_STAT_HEIGHT]),
                    im_width, im_height, FRAME_EDGE_MARGIN_PX,
                ))
                marks.append((int(centroids[i][0]), int(centroids[i][1]), int(np.sqrt(area_px / np.pi) * 1.1)))
                labels.append(i)
        return areas_um2, corrections, marks, labels

    def _particle_outlines(self, labels: np.ndarray, stats: np.ndarray, num_labels: int,
                           counted: Sequence[int]) -> Dict[str, list]:
        """Outlines of the particles of a frame (see particle_outlines): those counted (labels ``counted``, in the
        order of the frame's particles), and the others of at least min_par_pixels pixels."""
        from autoemx.core.particle_outlines import frame_outlines

        counted_set = set(counted)
        others = [i for i in range(1, num_labels) if i not in counted_set
                  and int(stats[i, cv2.CC_STAT_AREA]) >= self.powder_meas_cfg.min_par_pixels]
        return frame_outlines(labels, counted, others)

    def _sample_output_dir(self) -> Optional[str]:
        """Folder of the sample (results_dir, or its parent if results_dir is the images folder)."""
        if self.results_dir is None:
            return None
        if os.path.basename(self.results_dir) == cnst.IMAGES_DIR:
            return os.path.dirname(self.results_dir)
        return self.results_dir

    def _save_particle_outlines(self, segmentation_id: int, frame_id: str, outlines: Optional[Dict[str, list]]):
        """Save the particle outlines of a frame (shown in the GUI); a failure only skips them."""
        sample_dir = self._sample_output_dir()
        if outlines is None or sample_dir is None:
            return
        from autoemx.core.particle_outlines import save_outlines

        try:
            save_outlines(sample_dir, segmentation_id, frame_id, outlines)
        except OSError as exc:
            logger.warning(f"Particle outlines of frame {frame_id} not saved: {exc}")

    def segment_frame(self, frame_image: np.ndarray, pixel_size_um: float, d_range_um: Tuple[float, float],
                      segmentation_id: int, with_outlines: bool = False):
        """
        Particles of a saved frame image valid for particle size statistics, with the current segmentation settings,
        among the diameters ``d_range_um`` searched in the frame (raised to the smallest particles of at least
        min_par_pixels pixels): the result of segmentation ``segmentation_id`` of the frame (FrameSegmentation), and
        with ``with_outlines`` the outlines of its particles (see particle_outlines), as a tuple.
        """
        self._im_height, self._im_width = frame_image.shape[:2]
        self.EM.pixel_size_um = pixel_size_um
        d_lo = max(float(d_range_um[0]), min_valid_diameter_um(pixel_size_um, self.powder_meas_cfg.min_par_pixels))
        d_hi = float(d_range_um[1])
        self.powder_meas_cfg.min_area_par = diameter_um_to_area_um2(min(d_lo, d_hi))
        self.powder_meas_cfg.max_area_par = diameter_um_to_area_um2(d_hi)
        par_mask, _ = self._get_particles_on_substrate_mask(cv2.GaussianBlur(frame_image, (5, 5), 0), save_image=False)
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)
        areas, corrections, marks, counted = self._valid_particles(frame_image, stats, centroids, num_labels)
        result = FrameSegmentation(segmentation_id=segmentation_id, d_range_um=(min(d_lo, d_hi), d_hi),
                                   particle_areas_um2=areas, particle_edge_corrections=corrections,
                                   particle_marks_px=marks)
        if with_outlines:
            return result, self._particle_outlines(labels, stats, num_labels, counted)
        return result

    def _current_frame_position(self) -> Optional[Tuple[float, float]]:
        """Stage position (mm) of the current frame, or None if not known (e.g. not read after a manual move)."""
        navigator = getattr(self.EM, "frame_navigator", None)
        if navigator is None or not getattr(navigator, "current_pos_is_measured", True):
            return None
        return getattr(navigator, "_current_pos", None)

    def _record_eds_search_frame(
        self,
        frame_image_shape: Tuple[int, ...],
        num_labels: int,
        stats: np.ndarray,
        image_file: Optional[str] = None,
        centroids: Optional[np.ndarray] = None,
        labels: Optional[np.ndarray] = None,
    ) -> None:
        """
        Record a frame searched for EDS particles (with the outlines of its particles if ``labels``, the label image,
        is given), with the area, edge correction and mark of each particle valid
        for particle size statistics (within the area limits, of at least min_par_pixels pixels, not touching the
        frame edge), whether it is then measured by EDS or not; stored as a frame of the sample ledger (see
        ``LedgerModule._store_eds_search_frames``).
        """
        im_height, im_width = frame_image_shape[:2]
        pixel_size_um = float(self.EM.pixel_size_um)
        if centroids is None:
            centroids = np.array([[stats[i, cv2.CC_STAT_LEFT] + stats[i, cv2.CC_STAT_WIDTH] / 2,
                                   stats[i, cv2.CC_STAT_TOP] + stats[i, cv2.CC_STAT_HEIGHT] / 2]
                                  for i in range(num_labels)])
        areas_um2, corrections, marks, counted = self._valid_particles(np.empty((im_height, im_width)), stats,
                                                                       centroids, num_labels)
        frame_width_um = im_width * pixel_size_um
        d_min_um = max(
            float(2.0 * np.sqrt(self.powder_meas_cfg.min_area_par / np.pi)),
            min_valid_diameter_um(pixel_size_um, self.powder_meas_cfg.min_par_pixels),
        )
        d_max_um = float(2.0 * np.sqrt(self.powder_meas_cfg.max_area_par / np.pi))
        if d_min_um >= d_max_um:
            return  # frames too coarse to measure particles of this size
        position = self._current_frame_position()
        label = self.EM.current_frame_label
        self.pending_eds_frames.append((
            (d_min_um, d_max_um, frame_width_um),
            ParticleStatsFrame(
                frame_id=f"EDS-{label}" if label is not None else "EDS",
                center_mm=Coordinate2D(x=float(position[0]), y=float(position[1])) if position is not None else None,
                frame_width_um=frame_width_um,
                pixel_size_um=pixel_size_um,
                area_um2=float(im_height * im_width * pixel_size_um**2),
                n_particles=len(areas_um2),
                particle_areas_um2=areas_um2,
                particle_edge_corrections=corrections,
                image_file=image_file,
                particle_marks_px=[(x, y, r, i) for i, (x, y, r) in enumerate(marks)],
            ),
            self._particle_outlines(labels, stats, num_labels, counted) if labels is not None else None,
        ))

    #%% Particle Segmentation Operations
    # =============================================================================
    def _get_particles_coordinates_in_frame(self, frame_image=None, pixel_size=None, results_dir=None):
        '''
        Detects and extracts center coordinates and widths of suitable particles in the current frame.
    
        If running at the microscope, moves to the next frame and collects the image. If an image and pixel size 
        are provided, it uses these instead. The function applies a mask to segment particles, filters them 
        by area, and determines their positions and frame sizes. It also saves a visualization image with 
        circles drawn around detected particles for development purposes.
    
        Parameters
        ----------
        frame_image : ndarray, optional
            Grayscale image of the current frame. If not provided and running at the EM, the image is acquired.
        pixel_size : float, optional
            Pixel size in micrometers. Required if not running at the EM.
        results_dir : str, optional
            Directory to save results. Required for saving images in development mode.
    
        Returns
        -------
        bool or None
            True if at least one suitable particle was found in the frame.
            False if no particles were found.
            None if there are no more frames available.
        
        Notes
        -----
        - The function saves a visualization image with detected particles for development.
        - All commented-out imshow and development lines are preserved for debugging.
        - Particle positions and frame widths are stored as class attributes for later use.
        
        Potential Improvements
        ---------------------
        - Exclude particles that are close to larger particles in the direction of the EDS (Energy Dispersive Spectroscopy) detector.
          Large particles in the path can absorb or scatter X-rays emitted from smaller particles, degrading the quality and accuracy 
          of the spectral signal for those particles.
        
          Suggested implementation:
          - For each detected particle, determine if there are larger particles located "upstream" (i.e., between the particle and 
            the EDS detector direction).
          - To do this, add a control '_is_particle_shadowed()' together with '_is_particle_area_ok()' to select valid particles.
        '''
        if frame_image is not None and pixel_size is not None:
            # Image given (e.g. a saved frame reanalysed): used even at the microscope
            self._im_height, self._im_width = frame_image.shape
            self.EM.pixel_size_um = pixel_size
            if results_dir:
                self.results_dir = results_dir
        elif EM_driver.is_microscope_connected():
            self._check_EM_controller_initialization()
            move_to_frame_success = self.EM.go_to_next_frame()
            if move_to_frame_success:
                # Collect image
                frame_image = EM_driver.get_image_data(self._im_width, self._im_height, 1)
            else:
                # No more frames are available
                return None
        else:
            raise ValueError('Function "_get_particles_coordinates_in_frame()" must be run at the microscope, or it needs to be passed both image and its pixel size')
        
        
        ### Select particles on substrate
        # Apply a Gaussian blur to remove single bright pixels (the image as acquired is saved below)
        acquired_image = frame_image
        frame_image = cv2.GaussianBlur(frame_image, (5, 5), 0)
        
        # Get mask of particles on substrate
        par_mask, _ = self._get_particles_on_substrate_mask(frame_image)
    
        # Find connected components
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)
        
        ### Filter out particles that are too small or too big, and store their centroid + size
        par_pos_pixels = []
        par_fw_pixels = []
        par_radius_pixels = []
        fw_margin_um = 5 # Pixel size in frame is a few um, so the particle will be shifted
        # fw_scale_factor = 1.8 # Multiplicative factor to obtain margins from particle
        for i in range(1, num_labels):  # Skip the background component (index 0)
            if self._is_particle_area_ok(stats[i, cv2.CC_STAT_AREA]):
                
                # Append particle centroid
                par_pos_pixels.append(centroids[i])
                
                # Calculate what would be the frame width in pixels in order to contain the particle fully in its width
                fw_width = stats[i, cv2.CC_STAT_WIDTH]
                # Same as above, but to contain the particle in its height
                fw_height = stats[i, cv2.CC_STAT_HEIGHT] / self._im_height * self._im_width
                
                # Select the largest frame width to make sure it contains the particle fully
                par_fw = max([fw_width, fw_height])
                # Append to list of frame_widths
                par_fw_pixels.append(par_fw)
                
                # Save particle radius, for proper circle size in saved image
                par_radius_pixels.append(max(fw_width, stats[i, cv2.CC_STAT_HEIGHT]) / 2)
                
        ### Visualize selected particles. Used for code development
        #     else:
        #         # Blacken pixels corresponding to excluded components
        #         par_mask[labels == i] = 0
        # cv2.imshow('Filtered mask', par_mask)
        
        
        # Save the frame image analysed (annotated with the identified particles only if annotate_particle_images),
        # and count the particles of the frame for particle size statistics
        particle_marks = [
            (int(center[0]), int(center[1]), int(rad), None) for rad, center in zip(par_radius_pixels, par_pos_pixels)
        ]
        image_file = self._save_frame_images(
            self._unique_image_stem(f"{self._sample_id}_frEDS-{self.EM.current_frame_label}"), acquired_image,
            particle_marks,
        )
        self._record_eds_search_frame(frame_image.shape, num_labels, stats, image_file, centroids, labels)
        
        
        # Return false if no particles were detected in the frame
        num_par = len(par_pos_pixels)
        if self.verbose:
            if num_par == 1:
                par_string = 'particle was'
            else:
                par_string = 'particles were'
            logger.info(f"✅ {num_par} {par_string} found in current frame")
        if num_par == 0:
            return False
        
        ### Convert dimensions from pixels to mm
        # Calculate absolute position of particle within the EM, in mm
        par_pos_abs_mm = self.EM.convert_pixel_pos_to_mm(par_pos_pixels)
        # Convert the frame width to mm
        par_fw_mm = (np.array(par_fw_pixels) * self.EM.pixel_size_um + 2 * fw_margin_um) * 10**-3
        
        # Store particle information
        self.par_pos_abs_mm = par_pos_abs_mm
        self.par_fw_mm = par_fw_mm
        self._num_par_in_frame = num_par
        
        # Initialise counter to keep track at how many particles within the frame have been analysed
        self._fr_par_cntr = 0
        
        # Returns True if at least 1 particle was found in the frame, otherwise returns False
        return True


    def _get_particles_on_substrate_mask(self, frame_image, save_image=False, persist_mask=True):
        '''
        Generates a binary mask of detected particles on the substrate from the input frame image.

        This function applies a brightness threshold to the input image to segment particles, 
        finds contours, and fills inner contours to ensure particles are fully masked. 
        Optionally, the mask image can be saved to disk. The function returns the mask and 
        the path where the mask image would be saved.

        Parameters
        ----------
        frame_image : ndarray
            The grayscale input image of the current frame.
        save_image : bool, optional
            If True, saves the binary mask image to disk (default: False).
        persist_mask : bool, optional
            If False, skip writing the mask image even in development mode. Used when a
            caller (e.g. per-particle masking) will save its own mask under a different name.

        Returns
        -------
        par_mask : ndarray
            Either a binary mask of detected particles, or a labels array, where the pixels
            of each particle are identified by a different integer (same as labels returned by
                                                                    cv2.ConnectedComponents)
        mask_img_path : str
            File path for where the mask image is (or would be) saved.
            
        Note
        ----
        The commented `cv2.imshow` line can be enabled for debugging visualization.

        '''
        if self.powder_meas_cfg.par_segmentation_model not in self.powder_meas_cfg.AVAILABLE_PAR_SEGMENTATION_MODELS:
            self.powder_meas_cfg.par_segmentation_model = "threshold_bright"
            warnings.warn(
                f"Chosen particle segmentation model {self.powder_meas_cfg.par_segmentation_model} not available."
                "Defaulting to 'threshold_bright'",
                UserWarning
            )
            
        if self.powder_meas_cfg.par_segmentation_model == "threshold_bright":
            # Apply the threshold to get a binary image
            _, par_mask = cv2.threshold(frame_image, self.powder_meas_cfg.par_brightness_thresh, 255, cv2.THRESH_BINARY)
            
            # Find all contours in the image and fill them
            contours, hierarchy = cv2.findContours(par_mask, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
            for i, contour in enumerate(contours):
                if hierarchy[0][i][3] != -1:  # If contour is inside another contour
                    cv2.drawContours(par_mask, [contour], 0, 255, -1)
                    
        elif self.powder_meas_cfg.par_segmentation_model in self.powder_meas_cfg.AVAILABLE_PAR_SEGMENTATION_MODELS:
            model_module = par_seg_models.require_model(self.powder_meas_cfg.par_segmentation_model)
            par_mask = model_module.segment_particles(frame_image, self.powder_meas_cfg, save_image, self.EM)

        else:
            raise ValueError(f"Unknown error with current particle segmentation model {self.powder_meas_cfg.par_segmentation_model}")
        
        # cv2.imshow('Segmented Particles Mask', par_mask)
        
        mask_img_path = os.path.join(self.results_dir, self._sample_id + f'_fr{self.EM.current_frame_label}' + '_mask.png')
        # Mask is always saved when collecting particles. Avoids double saving
        save_image = save_image and not self.EM.measurement_cfg.type == self.EM.measurement_cfg.PARTICLE_STATS_MEAS_TYPE_KEY
        if persist_mask and (self.development_mode or save_image):
            # Scale bar drawn on a copy: the mask itself must stay binary for the particle detection. Label images
            # are saved as binary masks (cv2 does not write int32 images)
            mask_to_save = par_mask.copy() if is_binary_mask(par_mask) else np.where(par_mask > 0, 255, 0).astype(np.uint8)
            if mask_to_save.dtype == np.bool_:
                mask_to_save = mask_to_save.astype(np.uint8) * 255
            draw_scalebar(mask_to_save, self.EM.pixel_size_um)
            cv2.imwrite(mask_img_path, mask_to_save)
        
        return par_mask, mask_img_path
    
    
    def _get_connected_components_with_stats(self, par_mask: np.ndarray):
        """
        Compute connected components with statistics.
    
        This function accepts either:
          1. A binary mask (0/255 or boolean), in which case it directly calls
             cv2.connectedComponentsWithStats.
          2. A pre-labeled image (integer labels, like output of cv2.connectedComponents),
             in which case stats and centroids are recomputed manually.
    
        Parameters
        ----------
        par_mask : np.ndarray
            Input binary mask (bool, or uint8 with values 0/255) or label image (any other integer image: an
            int32 image with labels 0 and 255 only is a label image).
    
        Returns
        -------
        Same as cv2.connectedComponentsWithStats
        
        num_labels : int
            Number of connected components (including background).
        labels : np.ndarray
            Labeled image of the same size as input.
        stats : np.ndarray
            Statistics for each label. Shape: (num_labels, 5).
            Columns indexed by:
                cv2.CC_STAT_LEFT   (x)
                cv2.CC_STAT_TOP    (y)
                cv2.CC_STAT_WIDTH  (width)
                cv2.CC_STAT_HEIGHT (height)
                cv2.CC_STAT_AREA   (area)
        centroids : np.ndarray
            Centroids of each component. Shape: (num_labels, 2).
        """
    
        # --- Case 1: Binary image (bool, or uint8 with values 0 and 255 only) ---
        if is_binary_mask(par_mask):
            num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(
                par_mask.astype(np.uint8), connectivity=8, ltype=cv2.CV_32S
            )

        else:
            # --- Case 2: Already a label image ---
            labels = par_mask.astype(np.int32, copy=False)
            num_labels, stats, centroids = label_image_stats(labels)

        return num_labels, labels, stats, centroids

    
    def is_particle_at_frame_edge(self, stats, i):
        '''
        Determines whether a particle is located at or near the edge of the image frame.
    
        A margin is used to account for detection tolerances. If any part of the particle's bounding box 
        is within `pixel_margin` pixels of the image border, it is considered to be at the edge.
    
        Parameters
        ----------
        stats : ndarray
            Connected components statistics array (as returned by OpenCV), where each row corresponds 
            to a detected particle and columns to bounding box info (LEFT, TOP, WIDTH, HEIGHT, etc.).
        i : int
            Index of the particle in the stats array.
    
        Returns
        -------
        bool
            True if the particle is at or near the edge of the image frame, False otherwise.
        '''
        # Apply margin to account for detection tolerances
        pixel_margin = FRAME_EDGE_MARGIN_PX
    
        # Check if the particle's bounding box touches or is near any image edge
        is_particle_at_frame_edge = any([
            stats[i, cv2.CC_STAT_LEFT] <= pixel_margin,
            stats[i, cv2.CC_STAT_TOP] <= pixel_margin,
            stats[i, cv2.CC_STAT_LEFT] + stats[i, cv2.CC_STAT_WIDTH] >= self._im_width - pixel_margin,
            stats[i, cv2.CC_STAT_TOP] + stats[i, cv2.CC_STAT_HEIGHT] >= self._im_height - pixel_margin
        ])
    
        return is_particle_at_frame_edge


    def _is_particle_area_ok(self, area):
        '''
        Checks whether the given particle area is within acceptable limits.
    
        The function converts the minimum and maximum allowed particle areas from µm² to pixel units 
        using the current pixel size. It then determines if the provided area is within this range.
        If in manual collection mode and the area is out of bounds, it prints a message to the user.
    
        Parameters
        ----------
        area : float
            The area of the particle in pixels.
    
        Returns
        -------
        bool
            True if the particle area is within the allowed range, False otherwise.
        '''
        # Calculate acceptable particle area dimensions with current pixel size
        min_area_pixels = self.powder_meas_cfg.min_area_par / self.EM.pixel_size_um**2  # e.g., 4um^2 for a 2um x 2um particle
        max_area_pixels = self.powder_meas_cfg.max_area_par / self.EM.pixel_size_um**2  # e.g., 100um^2 for a 10um x 10um particle
    
        is_par_large_enough = area >= min_area_pixels
        if not is_par_large_enough and self.is_manual_particle_selection:
            logger.warning(f"⚠️ Selected particle is too small ({area*self.EM.pixel_size_um**2:.2f} um^2), select another one")
    
        is_par_small_enough = area <= max_area_pixels
        if not is_par_small_enough and self.is_manual_particle_selection:
            logger.warning(f"⚠️ Selected particle is too large ({area*self.EM.pixel_size_um**2:.2f} um^2), select another one")
    
        is_par_area_ok = is_par_large_enough and is_par_small_enough
    
        return is_par_area_ok


    def _capture_ref_frame(self, par_image=None, pixel_size_um=None):
        '''
        Capture (or accept) the grayscale frame used for spot selection.

        At the microscope, a fresh frame is acquired via the driver. Offline, the
        provided ``par_image`` is used and image dimensions / pixel size are
        updated from it. This is the single place that defines how the reference
        frame is obtained, shared by automated and callback spot selection.

        Parameters
        ----------
        par_image : ndarray, optional
            Grayscale frame to use when not running at the microscope.
        pixel_size_um : float, optional
            Pixel size in micrometers. Required together with ``par_image`` offline.

        Returns
        -------
        ndarray
            The grayscale frame image.
        '''
        if EM_driver.is_microscope_connected():
            self._check_EM_controller_initialization()
            return EM_driver.get_image_data(self._im_width, self._im_height, 1)
        if par_image is not None and pixel_size_um is not None:
            self._im_height, self._im_width = par_image.shape
            self.EM.pixel_size_um = pixel_size_um
            return par_image
        raise ValueError('This function must be run at the microscope, or it needs to be passed both image and its pixel size')


    def _segment_center_particle(
        self,
        par_image: np.ndarray,
        *,
        centering: bool = False,
    ) -> Optional[Tuple[np.ndarray, int, np.ndarray, np.ndarray, np.ndarray]]:
        """
        Segment the particle at the image center with the configured segmentation model.

        Returns
        -------
        tuple or None
            ``(binary_mask, par_label, stats, centroids, source_image)`` for the
            selected particle. ``source_image`` is the frame the mask was computed
            on and may differ from ``par_image`` if the stage was recentered.
            Area is computed from this full (non-eroded) mask before any spot-selection
            erosion. Returns None if no suitable particle is found.
        """
        par_mask, _ = self._get_particles_on_substrate_mask(par_image, persist_mask=False)
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)

        if num_labels == 1:
            return None

        center_x = int(self._im_width / 2)
        center_y = int(self._im_height / 2)
        par_label = int(labels[center_y, center_x])

        center_on_particle = par_label > 0
        center_area_ok = center_on_particle and self._is_particle_area_ok(
            stats[par_label, cv2.CC_STAT_AREA]
        )

        # Prefer the particle under the image center. Offline, keep it even when it
        # is outside configured area limits so size/spot selection still target the
        # particle being analysed. At the microscope, oversized/undersized center
        # particles trigger a re-center onto the nearest area-ok neighbour.
        if center_area_ok or (center_on_particle and not EM_driver.is_microscope_connected()):
            binary_mask = np.where(labels == par_label, 255, 0).astype(np.uint8)
            return binary_mask, par_label, stats, centroids, par_image

        if centering:
            return None

        distances = np.linalg.norm(centroids[1:] - np.array([center_x, center_y]), axis=1)
        sorted_indices = np.argsort(distances) + 1  # Skip background (index 0)
        for label in sorted_indices:
            label = int(label)
            if not EM_driver.is_microscope_connected():
                binary_mask = np.where(labels == label, 255, 0).astype(np.uint8)
                return binary_mask, label, stats, centroids, par_image
            if self._is_particle_area_ok(stats[label, cv2.CC_STAT_AREA]):
                new_center = self.EM.convert_pixel_pos_to_mm(centroids[label])
                self.EM.move_to_pos(new_center)
                center_xy = np.asarray(new_center, dtype=float).reshape(-1)
                self.current_par_coordinates_mm = (float(center_xy[0]), float(center_xy[1]))
                recentered_image = self._capture_ref_frame()
                return self._segment_center_particle(recentered_image, centering=True)
        return None

    def _store_current_particle_area(self, area_pixels: float) -> float:
        """Store particle area in µm² from a non-eroded pixel area."""
        area_um = float(area_pixels) * float(self.EM.pixel_size_um) ** 2
        self.current_par_area_um2 = area_um
        return area_um

    def _get_particle_mask(self, par_image=None, pixel_size_um=None, results_dir=None, centering=False):
        '''
        Returns a binary mask of the particle at the center of the image, if the particle is large enough.
        If no particle is present at the center, the function searches for and centers on the next closest 
        sufficiently large particle. This approach makes the function robust against drift or movement of 
        the particle during analysis. If no valid particle is found, the function returns None.
    
        Parameters
        ----------
        par_image : ndarray, optional
            Grayscale image of the current frame. If not provided and running at the EM, the image is acquired.
        pixel_size_um : float, optional
            Pixel size in micrometers. Required if not running at the EM.
        results_dir : str, optional
            Directory to save results. Required for saving the mask in development mode.
        centering : bool, optional
            If True, indicates this is a second attempt at centering on a particle (prevents infinite recursion).
    
        Returns
        -------
        tuple or None
            (par_mask, par_image): Binary mask of the selected particle and the corresponding image.
            None: If no suitable particle is found.
        
        Notes
        -----
        - Uses ``powder_meas_cfg.par_segmentation_model`` (default ``threshold_bright``).
        - Stores ``current_par_area_um2`` from the full mask before any edge erosion.
        - The function is robust to small misalignments: if the center particle is not valid, it finds the next closest.
        - The returned ``par_image`` is always the frame the mask was computed on, including after a microscope recenter.
        - In development mode, it saves the resulting mask with a scalebar overlay.
        - The commented `cv2.imshow` lines can be enabled for debugging visualization.
        '''
        if not EM_driver.is_microscope_connected() and results_dir:
            self.results_dir = results_dir
        par_image = self._capture_ref_frame(par_image, pixel_size_um)

        segmented = self._segment_center_particle(par_image, centering=centering)
        if segmented is None:
            if not centering:
                self.current_par_area_um2 = None
            return None

        par_mask, par_label, stats, _centroids, par_image = segmented
        # Area from the full (non-eroded) mask — spot selection may erode later.
        self._store_current_particle_area(stats[par_label, cv2.CC_STAT_AREA])

        if self.development_mode and self.results_dir:
            mask_to_save = draw_scalebar(par_mask.copy(), self.EM.pixel_size_um)
            cv2.imwrite(os.path.join(
                self.results_dir,
                self._sample_id + f'_par{self.tot_par_cntr}_fr{self.EM.current_frame_label}_mask.png'
            ), mask_to_save)

        return (par_mask, par_image)


    def _erode_particle_mask(self, par_mask, margin, erode_inner=False):
        """
        Erode a binary particle mask by a specified margin.
    
        This function erodes the mask either only at the outer boundary or at both the outer boundary and inner holes,
        depending on the `erode_inner` flag.
    
        Parameters
        ----------
        par_mask : np.ndarray
            Binary mask of the particle (dtype uint8, values 0 and 255).
        margin : int
            The erosion margin (in pixels). The structuring element will be of size (2 * margin + 1).
        erode_inner : bool, optional
            If False (default), only the outer boundary of the mask is eroded,
            leaving holes inside the particle unaffected.
            If True, both the outer boundary and any internal holes are eroded.
    
        Returns
        -------
        final_mask : np.ndarray
            The eroded binary mask (same dtype and shape as `par_mask`).
    
        Notes
        -----
        - This method assumes the mask is a binary 8-bit image (0 for background, 255 for foreground).
        - The erosion uses an elliptical structuring element.
    
        Examples
        --------
        >>> eroded = self._erode_particle_mask(par_mask, margin=3, erode_inner=False)
        """
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * margin + 1, 2 * margin + 1))
        if erode_inner:
            # Erode everywhere (outer and inner edges)
            final_mask = cv2.erode(par_mask, kernel)
        else:
            # Erode only the outer contour, not holes
            contours, hierarchy = cv2.findContours(par_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            outer_contour_mask = np.zeros_like(par_mask)
            cv2.drawContours(outer_contour_mask, contours, -1, 255, -1)  # Fill outer boundaries
            eroded_outer_contour_mask = cv2.erode(outer_contour_mask, kernel)
            final_mask = cv2.bitwise_and(par_mask, eroded_outer_contour_mask)
        return final_mask
    
    
    def _find_particle_bright_regions(self, final_mask, par_image):
        """
        Detect bright regions within a particle mask in the input image, in order to acquire X-ray spectra
        from thickest regions of particle and maximize X-ray intensity.
    
        This function applies a mask to the input image, blurs it, normalizes the intensity,
        thresholds to find bright spots, and returns the thresholded image and the minimum area (in pixels)
        for spot detection.
    
        Parameters
        ----------
        final_mask : np.ndarray
            Binary mask (uint8, values 0 and 255) specifying the region of interest (the particle).
        par_image : np.ndarray
            Grayscale image (uint8) of the particle region.
    
        Returns
        -------
        thresholded_image : np.ndarray
            Binary image (uint8, values 0 and 255) showing detected bright spots within the particle mask.
        min_area_pixels : int
            Minimum spot area in pixels, calculated as corresponding to 0.1 μm².
    
        Notes
        -----
        - The function applies a Gaussian blur before normalization and thresholding.
        - The threshold value is taken from `self.powder_meas_cfg.par_xy_spots_thresh`.
        - The minimum area is computed from `self.EM.pixel_size_um`.
        - If `self.development_mode` is True and `self.results_dir` is set, a debug image is saved with a scalebar.
    
        Examples
        --------
        >>> spots_mask, min_area = self._find_particle_bright_regions(final_mask, par_image)
        """
        masked_image = cv2.bitwise_and(par_image, par_image, mask=final_mask)
        # cv2.imshow('CV Image', masked_image)
        if np.all(masked_image == 0):
            return np.zeros_like(masked_image), 0
        masked_im_blurred = cv2.GaussianBlur(masked_image, (5, 5), 0)
        # cv2.imshow('Blurred Image', masked_im_blurred)
        
        # Normalise particle intensity to select thickest regions, regardless of intensity of neighbouring features
        norm_masked_image = (masked_im_blurred / np.max(masked_im_blurred) * 255).astype(np.uint8)
        _, thresholded_image = cv2.threshold(norm_masked_image, self.powder_meas_cfg.par_xy_spots_thresh, 255, cv2.THRESH_BINARY)
        # cv2.imshow('CV Image', thresholded_image)
        # cv2.imshow('CV Image 2', cv2.bitwise_and(par_image, thresholded_image))
        min_area_pixels = int(0.1 / self.EM.pixel_size_um ** 2)
        if self.development_mode and self.results_dir:
            eroded_par_mask = draw_scalebar(thresholded_image.copy(), self.EM.pixel_size_um)
            cv2.imwrite(os.path.join(self.results_dir, self._sample_id + f'_par{self.tot_par_cntr}_fr{self.EM.current_frame_label}_maskXY.png'), eroded_par_mask)
        return thresholded_image, min_area_pixels


    def prepare_mask_for_visualization(self, mask: np.ndarray) -> np.ndarray:
        """
        Prepare a segmentation mask for visualization.
    
        Behavior:

        - Binary masks with values in {0, 1} are scaled to [0, 255].
        - Binary masks with values in {0, 255} are returned unchanged.
        - Integer label masks have brightness reversed so higher labels are brighter,
          with a minimum brightness of 30 for positive values. Background (0) remains black.
        - Floating-point masks are normalized to [0, 255] with the same rules for
          positive values and background.
        """
        mask = np.asarray(mask)
        unique_vals = np.unique(mask)
    
        # --- Binary mask [0, 1] ---
        if np.array_equal(unique_vals, [0, 1]):
            return (mask * 255).astype(np.uint8)
    
        # --- Binary mask [0, 255] ---
        elif np.array_equal(unique_vals, [0, 255]):
            return mask.astype(np.uint8)
    
        # --- Integer label masks ---
        elif np.issubdtype(mask.dtype, np.integer):
            max_val = mask.max()
            if max_val > 0:
                scaled = mask.astype(np.float32)
                pos_mask = scaled > 0
    
                # Reverse intensity so higher labels → brighter
                scaled[pos_mask] = (max_val - scaled[pos_mask]) * (255.0 / max_val)
    
                # Minimum brightness for positive labels
                scaled[pos_mask] = np.clip(scaled[pos_mask], 30, 255)
                scaled[~pos_mask] = 0
                return scaled.astype(np.uint8)
            else:
                return mask.astype(np.uint8)
    
        # --- Floating-point masks ---
        else:
            min_val = mask.min()
            max_val = mask.max()
            if max_val > min_val:
                scaled = (mask - min_val) / (max_val - min_val) * 255.0
                scaled = scaled.astype(np.float32)
    
                pos_mask = mask > 0
                scaled[pos_mask] = np.clip(scaled[pos_mask], 30, 255)
                scaled[~pos_mask] = 0
    
                return scaled.astype(np.uint8)
            else:
                return np.zeros_like(mask, dtype=np.uint8)
    
    
    #%% Selection of spots for X-Ray spectra acquisition
    # =============================================================================
    def get_XS_acquisition_spots_coord_list(
            self, n_tot_sp_collected, 
            par_image=None, pixel_size_um=None, results_dir=None):
        '''
        Returns a list of pixel coordinates (in the current image) for X-ray spectrum spot collection on a particle.
    
        The function finds a suitable particle mask, erodes the mask to avoid particle edges, finds bright regions (or peak spots),
        and selects up to `powder_meas_cfg.max_spectra_per_par` spots per particle with a minimum distance between them (determined through 'powder_meas_cfg.par_mask_margin').
        The selection strategy for features and spacing can be controlled via powder_meas_cfg.par_feature_selection and powder_meas_cfg.par_spot_spacing.
    
        Parameters
        ----------
        n_tot_sp_collected : int
            Counter for the total number of spectra collected (used for labeling).
        par_image : ndarray, optional
            Grayscale image of the current frame. If not provided and running at the SEM, the image is acquired.
        pixel_size_um : float, optional
            Pixel size in micrometers. Required if not running at the SEM.
        results_dir : str, optional
            Directory to save result images.
    
        Returns
        -------
        selected_points : list[tuple[int, int]]
            List of selected (x, y) spot coordinates in image pixel units.
    
        Potential Improvements
        ---------------------
        - Select points to ensure all phases are individuated, biasing spot selection to ensure representation
            of phases with different bright/dark contrast. Currently, only the brightest spots are picked,
            which may miss some phases.
        - Add detection of both peaks and valleys to ensure both bright and dark spots are tested.
        - Exclude points that are close to larger particles along the X-ray emission path to the EDS detector,
          as these may absorb emitted X-rays and degrade spectral signal.
        '''
        # In callback mode, spots are supplied by xsp_spot_selector. Skip the
        # particle mask computation entirely and reuse a single reference frame
        # per particle (see _get_callback_xsp_spots).
        if self.powder_meas_cfg.par_spot_selection_mode == 'callback':
            return self._get_callback_xsp_spots(
                n_tot_sp_collected,
                par_image=par_image,
                pixel_size_um=pixel_size_um,
            )

        # --- 1. Acquire or prepare the particle mask and image ---
        if EM_driver.is_microscope_connected():
            par_mask_return = self._get_particle_mask()
        elif par_image is not None and pixel_size_um is not None:
            self._im_height, self._im_width = par_image.shape
            self.EM.pixel_size_um = pixel_size_um
            par_mask_return = self._get_particle_mask(par_image, pixel_size_um)
            if results_dir:
                self.results_dir = results_dir
        else:
            raise ValueError('This function must be run at the microscope, or it needs to be passed both image and its pixel size')
    
        # Check if a particle was detected. If not, return empty list
        if par_mask_return is None:
            self.ref_image = None
            return []
        else:
            par_mask, par_image = par_mask_return
    
        # --- 2. Erode the particle mask to avoid edge effects ---
        margin = max(10, int(self.powder_meas_cfg.par_mask_margin / self.EM.pixel_size_um))
        final_mask = self._erode_particle_mask(par_mask, margin)
    
        # --- 3. Find bright points in image, which indicate highest regions on particle---
        thresholded_image, min_area_pixels = self._find_particle_bright_regions(final_mask, par_image)
    
        # --- 4. Collect candidate points from thresholded image ---
        all_points = self._collect_candidate_points(
            thresholded_image, par_image, self.powder_meas_cfg.par_feature_selection, min_area_pixels
        )
    
        if len(all_points) == 0:
            self.ref_image = None
            return []
    
        # --- 5. Select points with minimum distance and maximum count ---
        min_distance_xsp_spots = (self.powder_meas_cfg.xsp_spots_distance_um / self.EM.pixel_size_um)
        if self.powder_meas_cfg.par_spot_spacing == 'random':
            selected_points = EM_Particle_Finder._select_XSspots_randomly(
                all_points, max_points=self.powder_meas_cfg.max_spectra_per_par, min_distance=min_distance_xsp_spots
            )
        elif self.powder_meas_cfg.par_spot_spacing == 'maximized':
            selected_points = EM_Particle_Finder._select_evenly_spaced_XSspots(
                all_points, max_points=self.powder_meas_cfg.max_spectra_per_par, min_distance=min_distance_xsp_spots
            )
    
        # --- 6. Annotate image and save ---
        if self.development_mode and self.results_dir:
            color_image = cv2.cvtColor(par_image, cv2.COLOR_GRAY2BGR)
            for center in selected_points:
                # Add circle indicating where X-ray spectrum was collected
                cv2.circle(color_image, center, 10, (0, 0, 255), -1)
                label_pos = (center[0] - 30, center[1] - 15)
                cv2.putText(color_image, str(n_tot_sp_collected), label_pos, cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2, cv2.LINE_AA)
                n_tot_sp_collected += 1
            color_image = draw_scalebar(color_image, self.EM.pixel_size_um)
            # cv2.imshow('Selected XS spots', color_image)
            cv2.imwrite(os.path.join(self.results_dir, self._sample_id + f'_par{self.tot_par_cntr}_fr{self.EM.current_frame_label}_xyspots.png'), color_image)

        self.ref_image = par_image
        return [tuple(map(int, pt)) for pt in selected_points]

    def _get_callback_xsp_spots(
        self,
        n_tot_sp_collected: int,
        par_image=None,
        pixel_size_um=None,
    ) -> List[tuple]:
        """Invoke xsp_spot_selector and validate returned pixel coordinates.

        The reference frame is captured once per particle (when ``self.ref_image``
        is None) and reused for every spot batch on that particle, so the
        drift-correction reference stays fixed across that particle's spectra. No
        particle mask is computed in this mode.
        """
        if self.ref_image is None:
            self.ref_image = self._capture_ref_frame(par_image, pixel_size_um)
            # Measure particle size independently of the callback. Uses the configured
            # segmentation model on the raw (non-eroded) reference frame. If the
            # stage was recentered, keep the recaptured frame as the reference.
            segmented = self._segment_center_particle(self.ref_image)
            if segmented is not None:
                _mask, par_label, stats, _centroids, self.ref_image = segmented
                self._store_current_particle_area(stats[par_label, cv2.CC_STAT_AREA])
            else:
                self.current_par_area_um2 = None

        frame_label = getattr(self.EM, 'current_frame_label', '')
        if frame_label is None:
            frame_label = ''

        context = XSpSpotSelectionContext(
            particle_id=self.tot_par_cntr,
            n_tot_sp_collected=n_tot_sp_collected,
            ref_image=self.ref_image,
            pixel_size_um=float(self.EM.pixel_size_um),
            frame_label=str(frame_label),
            im_width=int(self._im_width),
            im_height=int(self._im_height),
            particle_finder=self,
        )
        try:
            raw_spots = self.xsp_spot_selector(context)
        except Exception as exc:
            logger.error(f"❌ xsp_spot_selector callback failed: {exc}")
            raise

        if not raw_spots:
            return []

        selected_points = validate_xsp_spot_pixels(
            raw_spots,
            im_width=int(self._im_width),
            im_height=int(self._im_height),
        )
        return selected_points

        
    def _collect_candidate_points(self, thresholded_image, par_image, feature_selection, min_area_pixels):
        """
        Collect candidate EDS/WDS spot coordinates from thresholded regions.
    
        This function identifies connected components (regions) in the thresholded image,
        filters them by minimum area, and collects candidate points for X-ray spot acquisition.
        The method of selection depends on the `feature_selection` argument:
        - 'random': all pixel coordinates within the component are collected.
        - 'peaks': only the brightest pixel within each component is selected.
    
        Parameters
        ----------
        thresholded_image : np.ndarray
            Binary image (uint8, values 0 and 255) indicating candidate regions for spot selection.
        par_image : np.ndarray
            Grayscale image (uint8) of the particle, used for peak detection.
        feature_selection : str
            'random' to select all pixels in the region, 'peaks' to select local maxima (bright spots).
        min_area_pixels : int
            Minimum area (in pixels) for a component to be considered.
    
        Returns
        -------
        all_points : list of tuple
            List of (x, y) coordinates of candidate points (in pixel units).
    
        Notes
        -----
        - For 'random', all pixels in each sufficiently large component are returned.
        - For 'peaks', local maxima are found by repeatedly masking out regions around each peak.
        - The margin for masking out peaks is determined by `self.powder_meas_cfg.par_mask_margin` and `self.EM.pixel_size_um`.
        - The maximum number of iterations for peak finding is `100 * self.powder_meas_cfg.max_spectra_per_par`.
        - Only components with area >= `min_area_pixels` are considered.
        """
        all_points = []
        num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(thresholded_image, connectivity=8)
        margin = max(10, int(self.powder_meas_cfg.par_mask_margin / self.EM.pixel_size_um))
        for i in range(1, num_labels):  # Skip background
            if stats[i, cv2.CC_STAT_AREA] >= min_area_pixels:
                if feature_selection == 'random':
                    components_coords = np.where(labels == i)
                    all_points += list(zip(components_coords[1], components_coords[0]))
                elif feature_selection == 'peaks':
                    component_mask = (labels == i).astype(np.uint8) * 255
                    component_image = cv2.bitwise_and(par_image, par_image, mask=component_mask)
                    component_image = cv2.GaussianBlur(component_image, (5, 5), 0)
                    max_iter = 100 * self.powder_meas_cfg.max_spectra_per_par
                    iter_cntr = 0
                    while iter_cntr < max_iter:
                        iter_cntr += 1
                        min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(component_image)
                        if max_val < self.powder_meas_cfg.par_xy_spots_thresh:
                            break
                        _, thresh = cv2.threshold(component_image, int(max_val * 0.97), 255, cv2.THRESH_BINARY)
                        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                        for contour in contours:
                            x, y, w, h = cv2.boundingRect(contour)
                            center_x = x + w // 2
                            center_y = y + h // 2
                            center_pnt = (center_x, center_y)
                            all_points.append(center_pnt)
                            cv2.circle(component_image, center_pnt, margin, (0, 0, 0), -1)
        return all_points
    
    
    @staticmethod
    def _select_XSspots_randomly(points, max_points, min_distance):
        '''
        Randomly selects up to `max_points` from a list of acquisition spots, ensuring that each selected point
        is at least `min_distance` away from all others already selected.
    
        Parameters
        ----------
        points : list of tuple or ndarray
            List of (x, y) coordinates to choose from.
        max_points : int
            Maximum number of points to select.
        min_distance : float
            Minimum allowed distance between any two selected points.
    
        Returns
        -------
        selected_points : list
            List of selected (x, y) points.
        '''
        # Shuffle points to unbias the choice of position
        random.shuffle(points)
    
        # Check if some selected points must be cut out
        if len(points) <= max_points:
            return points
    
        selected_points = [points[0]]
        for point in points:
            # Check if enough points have been selected
            if len(selected_points) >= max_points:
                break
    
            # Check if point is distant enough from all selected points
            if point not in selected_points and all(np.linalg.norm(np.array(point) - np.array(p)) > min_distance for p in selected_points):
                selected_points.append(point)
    
        return selected_points
    
    
    @staticmethod
    def _select_evenly_spaced_XSspots(points, max_points, min_distance):
        '''
        Selects up to `max_points` from a list of acquisition spots, maximizing the minimum distance between any two selected points.
        Points are chosen iteratively: each new point is the one farthest (in minimum distance) from those already selected.
    
        Parameters
        ----------
        points : list of tuple or ndarray
            List of (x, y) coordinates to choose from.
        max_points : int
            Maximum number of points to select.
        min_distance : float
            Minimum allowed distance between any two selected points (not strictly enforced, but likely achieved).
    
        Returns
        -------
        selected_points : list
            List of selected (x, y) points.
        '''
        # Shuffle points to unbias the choice of position
        random.shuffle(points)
    
        # Check if some selected points must be cut out
        if len(points) <= max_points:
            return points
    
        selected_points = [points[0]]
        while len(selected_points) < max_points:
            max_min_distance = -1                   
            best_point = None
            for point in points:
                if point not in selected_points:
                    min_dist_to_selected = min(np.linalg.norm(np.array(point) - np.array(p)) for p in selected_points)
                    if min_dist_to_selected > max_min_distance:
                        max_min_distance = min_dist_to_selected
                        best_point = point
            if best_point is not None:
                selected_points.append(best_point)
            else:
                break
        return selected_points
    
    #%% Particle Statistics
    # =============================================================================
    def attach_stats_ledger(self, ledger, ledger_path: Optional[str] = None) -> None:
        """Attach a SampleLedger whose ``particles`` list stores PSD ParticleInfo records."""
        self._ledger = ledger
        self._ledger_path = ledger_path
        if ledger is not None:
            self.analyzed_pars = ledger.particles

    def _persist_stats_ledger(self) -> None:
        """Write the attached particle-stats ledger to disk when available."""
        if self._ledger is None:
            return
        ledger_path = self._ledger_path
        if ledger_path is None and self.results_dir is not None:
            output_dir = (
                os.path.dirname(self.results_dir)
                if os.path.basename(self.results_dir) == cnst.IMAGES_DIR
                else self.results_dir
            )
            ledger_path = os.path.join(output_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
        if ledger_path is not None:
            self._ledger.to_json_file(ledger_path)

    def _set_size_band_area_limits(self, d_min_um: float, d_max_um: float, frame_width_um: float) -> None:
        """Set the particle area limits (and search frame width) of one size band."""
        self.powder_meas_cfg.min_area_par = diameter_um_to_area_um2(d_min_um)
        self.powder_meas_cfg.max_area_par = diameter_um_to_area_um2(d_max_um)
        self.powder_meas_cfg.par_search_frame_width_um = frame_width_um

    def _persist_stats_progress(self) -> None:
        """Save the particle size statistics measurement in progress (its frames) to the attached ledger."""
        if self._ledger is not None and self._stats_progress is not None:
            self._persist_stats_ledger()

    def _save_frame_images(
        self,
        filename: str,
        frame_image: Optional[np.ndarray],
        particle_marks: List[Tuple[int, int, int, Optional[int]]],
        mask_image: Optional[np.ndarray] = None,
    ) -> str:
        """
        Save a frame image (``<filename>_particles``) and optionally its particle mask (``<filename>_par_mask``),
        with the particles marked and a scale bar only if ``annotate_particle_images`` (they can be annotated later
        from the ledger, see ``autoemx.runners.annotate_particle_images``). Returns the frame image file name.
        """
        annotate = bool(getattr(self.EM.measurement_cfg, "annotate_particle_images", False))
        annotations = particle_mark_annotations(particle_marks, self._im_width) if annotate else None
        self.EM.save_frame_image(filename + '_particles', im_annotations=annotations, scalebar=annotate,
                                 frame_image=frame_image)
        if mask_image is not None:
            self.EM.save_frame_image(filename + '_par_mask', im_annotations=annotations, scalebar=annotate,
                                     frame_image=mask_image)
        return f"{filename}_particles.{self.EM.measurement_cfg.saved_images_extension}"

    def _record_stats_frame(
        self,
        frame_image_shape: Tuple[int, ...],
        frame_id: str,
        areas_um2: List[float],
        corrections: List[float],
        marks: List[Tuple[int, int, int]],
        image_file: Optional[str] = None,
        outlines: Optional[Dict[str, list]] = None,
    ) -> None:
        """
        Record the frame just analysed, with its particles, in the size band being collected and as a frame of the
        attached ledger (with the result of the segmentation of the measurement, and the outlines of its particles),
        saving the progress regularly.
        """
        band = self._stats_band
        if band is None:
            return
        im_height, im_width = frame_image_shape[:2]
        pixel_size_um = float(self.EM.pixel_size_um)
        position = self._current_frame_position()
        center = Coordinate2D(x=float(position[0]), y=float(position[1])) if position is not None else None
        frame = ParticleStatsFrame(
            frame_id=frame_id, center_mm=center, frame_width_um=im_width * pixel_size_um, pixel_size_um=pixel_size_um,
            area_um2=float(im_height * im_width * pixel_size_um**2), n_particles=len(areas_um2),
            particle_areas_um2=list(areas_um2), particle_edge_corrections=list(corrections), image_file=image_file,
            particle_marks_px=[(x, y, r, i) for i, (x, y, r) in enumerate(marks)],
        )
        band.add_frame(frame)
        progress = self._stats_progress
        if self._ledger is not None and progress is not None and progress.measurement is not None:
            window = (band.d_min_um, band.d_max_um)
            record = self._ledger.add_frame(FrameRecord(
                frame_id=frame_id, label=str(self.EM.current_frame_label), purpose="size_stats",
                measurement_id=progress.measurement.measurement_id, center_mm=center,
                frame_width_um=frame.frame_width_um, pixel_size_um=pixel_size_um, area_um2=frame.area_um2,
                search_window_um=window, image_file=image_file,
                segmentations=[FrameSegmentation(segmentation_id=progress.segmentation_id, d_range_um=window,
                                                 particle_areas_um2=list(areas_um2),
                                                 particle_edge_corrections=list(corrections),
                                                 particle_marks_px=list(marks))],
            ))
            progress.measurement.bands[progress.bands.index(band)].frame_ids.append(record.frame_id)
            self._save_particle_outlines(progress.segmentation_id, record.frame_id, outlines)
        if band.n_frames % STATS_SAVE_EVERY_N_FRAMES == 0:
            self._persist_stats_progress()

    def _add_stats_band(self, progress: StatsProgress, band: ParticleSizeBand) -> None:
        """Start a size band of the measurement (also in its ledger record)."""
        progress.bands.append(band)
        if progress.measurement is not None:
            progress.measurement.bands.append(SizeBand(d_min_um=band.d_min_um, d_max_um=band.d_max_um,
                                                       frame_width_um=band.frame_width_um))

    def _new_stats_progress(self, powder_meas_cfg: PowderMeasurementConfig, n_par_target: int,
                            count_eds_frames: bool = False) -> StatsProgress:
        """A new measurement, recorded in the attached ledger with the segmentation of its settings."""
        if self._ledger is None:
            return StatsProgress(n_par_target, powder_meas_cfg, [])
        segmentation_id = self._ledger.append_segmentation(
            SegmentationConfig.from_powder_meas_cfg(powder_meas_cfg, 0))
        measurement = SizeMeasurement(
            measurement_id=self._ledger.next_measurement_id(), n_par_target=n_par_target,
            powder_meas_cfg=powder_meas_cfg, segmentation_id=segmentation_id, count_eds_frames=count_eds_frames,
        )
        self._ledger.size_measurements.append(measurement)
        return StatsProgress(n_par_target, powder_meas_cfg, [], measurement, segmentation_id)

    def _collect_size_band(self, band: ParticleSizeBand, n_par_target: int) -> None:
        """Analyse frames until ``band`` holds ``n_par_target`` particles or no frame is left, recording each frame."""
        self._stats_band = band
        try:
            while band.n_particles < n_par_target:
                if self._stats_max_n_frames is not None and self._stats_progress is not None and sum(
                    b.n_frames for b in self._stats_progress.bands
                ) >= self._stats_max_n_frames:
                    break
                n_par_before = band.n_particles
                par_were_found = self._move_and_get_particles_stats_in_frame()
                if par_were_found is None:
                    logger.warning(
                        f"⚠️ Could not find {n_par_target} particles in this size band. "
                        f"Completed it with {band.n_particles} particles."
                    )
                    break
                if par_were_found is False:
                    if self.verbose:
                        logger.warning("⚠️ No particle was found in this frame")
                elif self.verbose:
                    n_par_found_frame = band.n_particles - n_par_before
                    n_par_band = band.n_particles
                    logger.info(
                        f"{n_par_found_frame} particle{'s' if n_par_found_frame != 1 else ''} "
                        f"{'were' if n_par_found_frame != 1 else 'was'} found in this frame."
                    )
                    logger.info(
                        f"A total of {n_par_band} particle{'s' if n_par_band != 1 else ''} "
                        f"{'have' if n_par_band != 1 else 'has'} now been analyzed in this size band. "
                        f"{max(n_par_target - n_par_band, 0)} more to go."
                    )
        finally:
            self._stats_band = None
        self._persist_stats_progress()

    def _resume_stats_progress(
        self, powder_meas_cfg: PowderMeasurementConfig, n_par_target: int
    ) -> Optional[StatsProgress]:
        """Measurement interrupted in the attached ledger, to resume if it has the same target and settings."""
        if self._ledger is None:
            return None
        measurement = next((m for m in reversed(self._ledger.size_measurements) if m.status == "in_progress"), None)
        if measurement is None:
            return None
        if measurement.n_par_target == n_par_target and measurement.powder_meas_cfg == powder_meas_cfg:
            bands = particle_statistics.size_stats_inputs(
                self._ledger, measurement.measurement_id, measurement.segmentation_id, None, False).bands
            for band, recorded in zip(bands, measurement.bands):  # the bands as planned by the measurement
                band.d_min_um = recorded.d_min_um
            progress = StatsProgress(n_par_target, powder_meas_cfg, bands, measurement, measurement.segmentation_id)
            n_frames = sum(band.n_frames for band in progress.bands)
            n_par = sum(band.n_particles for band in progress.bands)
            logger.info(
                f"ℹ️ Resuming the interrupted particle size statistics measurement "
                f"({n_frames} frames and {n_par} particles already analysed)."
            )
            return progress
        logger.warning(
            "⚠️ An interrupted particle size statistics measurement with different settings is not resumed: "
            "a new one starts (the frames of the interrupted one stay in the ledger)."
        )
        measurement.status = "interrupted"
        return None

    def _frame_width_range_um(self) -> Optional[Tuple[float, float]]:
        """Smallest and largest frame width of the microscope, in µm (None if not known)."""
        try:
            if EM_driver.is_microscope_connected():
                min_fw_mm, max_fw_mm = EM_driver.get_range_frame_width()
                return float(min_fw_mm) * 1000.0, float(max_fw_mm) * 1000.0
        except Exception:
            pass
        return None

    def _plan_band_frames(
        self, navigator, progress: StatsProgress, band_idx: int, fw_mm: float, reuse_grid: bool = False,
        d_range_um: Optional[Tuple[float, float]] = None,
    ) -> bool:
        """
        Plan the frames of a size band (diameters ``d_range_um``): the sample grid (first band) or subframes of the
        frames recorded in the previous band. When counting the frames searched by EDS (see
        :meth:`get_particle_stats`), these are also parent frames of the second band, and the frames overlapping
        an EDS frame searching the same sizes are left out. Returns False if the previous band has no frame.
        """
        if band_idx == 0:
            if (
                reuse_grid
                and navigator.grid_search_fw_mm is not None
                and np.isclose(float(navigator.grid_search_fw_mm), fw_mm, rtol=1e-6, atol=1e-9)
            ):
                navigator.reset_frame_cursor()
            else:
                navigator.configure_top_level_particle_frames(fw_mm)
        else:
            parent_frames = self._band_frame_descriptors(progress.bands[band_idx - 1])
            if band_idx == 1:
                # Frames searched by EDS: searched for the sizes of the first band, not for smaller ones
                parent_frames += [
                    (f.frame_id, (f.center_um[0] / 1000, f.center_um[1] / 1000), f.frame_width_um / 1000)
                    for f in self._eds_frames
                ]
            if not parent_frames:
                logger.warning("⚠️ No parent frames from the previous size band; skipping remaining finer bands.")
                return False
            navigator.configure_particle_subframes(parent_frames, fw_mm)
        if self._eds_frames and d_range_um is not None:
            self._skip_eds_searched_frames(navigator, d_range_um)
        return True

    def _skip_eds_searched_frames(self, navigator, d_range_um: Tuple[float, float]) -> None:
        """Leave out the planned frames overlapping a frame searched by EDS for the same particle sizes."""
        searched = [f for f in self._eds_frames if particle_statistics.size_ranges_overlap([f.d_range_um], d_range_um)[0]]
        n_skipped = navigator.skip_frames_overlapping([
            ((f.center_um[0] / 1000, f.center_um[1] / 1000), (f.half_size_um[0] / 1000, f.half_size_um[1] / 1000))
            for f in searched
        ])
        if n_skipped and self.verbose:
            logger.info(f"ℹ️ {n_skipped} frame{'s' if n_skipped != 1 else ''} already searched by EDS left out.")

    def _band_target(self, band: ParticleSizeBand, band_idx: int, n_par_target: int) -> int:
        """Particles still to find in a size band: the target, less the particles counted in the frames searched by
        EDS (when counted, see :meth:`get_particle_stats`) with a diameter in the band."""
        d = self._eds_diameters_um
        if not d.size:
            return n_par_target
        in_band = (d >= band.d_min_um * (1 - 1e-9)) & (
            d <= band.d_max_um * (1 + 1e-9) if band_idx == 0 else d < band.d_max_um * (1 - 1e-9)
        )
        n_eds = int(np.sum(in_band))
        if self.verbose:
            logger.info(f"ℹ️ {n_eds} particle{'s' if n_eds != 1 else ''} of this size band already counted in the "
                        f"frames searched by EDS: {max(n_par_target - n_eds, 0)} more to find.")
        return max(n_par_target - n_eds, 0)

    def _collect_in_band(self, band: ParticleSizeBand, band_idx: int, n_par_target: int) -> None:
        """Set the particle area limits of a size band and collect its particles."""
        self._set_size_band_area_limits(band.d_min_um, band.d_max_um, band.frame_width_um)
        if self.verbose:
            print_double_separator(
                f"Size band {band_idx + 1}: {band.d_min_um:.4g}–{band.d_max_um:.4g} µm (FOV {band.frame_width_um:.3g} µm)"
            )
        self._collect_size_band(band, n_par_target)

    def _continue_size_band(
        self, navigator, progress: StatsProgress, band_idx: int, n_par_target: int
    ) -> bool:
        """Collect more particles in a recorded size band, in its frames not analysed yet. False if none planned."""
        band = progress.bands[band_idx]
        if not self._plan_band_frames(navigator, progress, band_idx, band.frame_width_um / 1000.0,
                                      d_range_um=(band.d_min_um, band.d_max_um)):
            return False
        # Frames analysed already, by their navigation label
        records = [self._ledger.frame(f.frame_id) if self._ledger is not None else None for f in band.frames]
        navigator.skip_frames([r.label if r is not None else f.frame_id for f, r in zip(band.frames, records)])
        self._collect_in_band(band, band_idx, n_par_target)
        return True

    def _stats_eds_segmentation(self, progress: StatsProgress) -> Optional[int]:
        """Segmentation of the frames searched by EDS counted by the measurement (None: not counted)."""
        if self._ledger is None or not progress.powder_meas_cfg.include_eds_auto_particles:
            return None
        return particle_statistics.eds_segmentation_for(self._ledger, progress.segmentation_id)

    def _stats_eds_scans(self, progress: StatsProgress) -> List[ParticleSizeBand]:
        """Frames searched by EDS on the sample (in the attached ledger) to include in the statistics."""
        segmentation_id = self._stats_eds_segmentation(progress)
        if segmentation_id is None:
            return []
        return particle_statistics.size_stats_inputs(self._ledger, eds_segmentation_id=segmentation_id).eds_scans

    def _stats_particles(self, progress: StatsProgress, eds_scans: List[ParticleSizeBand]) -> List[ParticleInfo]:
        """Particles of the statistics besides those of frames: measured by EDS, without frames searched by EDS."""
        has_eds_frames = self._ledger is not None and bool(self._ledger.frames_of("eds_search"))
        return particle_statistics.select_size_stats_particles(
            self.analyzed_pars, progress.powder_meas_cfg.include_eds_auto_particles and not has_eds_frames,
        )

    def _size_stats_precision(self, progress: StatsProgress) -> Dict[str, float]:
        """Relative half-width of the 95 % confidence interval of the statistics whose precision is checked."""
        cfg = progress.powder_meas_cfg
        eds_scans = self._stats_eds_scans(progress)
        particles = self._stats_particles(progress, eds_scans)
        obs, weights, _, _ = particle_statistics.weighted_size_observations(particles, progress.bands, (), eds_scans)
        if not np.any(weights > 0):
            return {name: float("inf") for name in cfg.size_stats_ci_quantities}
        return particle_statistics.ci_relative_half_widths(
            particle_statistics.compute_particle_size_stats(obs.eq_diameters_um, weights),
            particle_statistics.bootstrap_confidence_intervals(particles, progress.bands, eds_scans=eds_scans),
            cfg.size_stats_ci_quantities,
        )

    def _collect_until_precise(self, progress: StatsProgress, navigator) -> bool:
        """
        Collect more particles in every size band (about 50 % more per round) until the confidence intervals of
        ``powder_meas_cfg.size_stats_ci_quantities`` are within ``size_stats_ci_tolerance``, the frame or time
        limit is reached, or no frame is left. Returns whether the precision was reached.
        """
        cfg = progress.powder_meas_cfg
        start_time = time.time()
        self._stats_max_n_frames = cfg.size_stats_max_n_frames
        try:
            while True:
                half_widths = self._size_stats_precision(progress)
                summary = ", ".join(f"{name} ±{100 * hw:.1f} %" for name, hw in half_widths.items())
                if all(hw <= cfg.size_stats_ci_tolerance for hw in half_widths.values()):
                    logger.info(f"✅ Precision reached: {summary} (95 % confidence).")
                    return True
                n_frames = sum(band.n_frames for band in progress.bands)
                elapsed_min = (time.time() - start_time) / 60
                if n_frames >= cfg.size_stats_max_n_frames or (
                    cfg.size_stats_max_time_min is not None and elapsed_min >= cfg.size_stats_max_time_min
                ):
                    logger.warning(f"⚠️ Frame or time limit reached before the precision: {summary}.")
                    return False
                if self.verbose:
                    logger.info(f"ℹ️ Precision {summary} (target ±{100 * cfg.size_stats_ci_tolerance:.1f} %): "
                                "collecting more particles.")
                n_frames_before = n_frames
                for band_idx, band in enumerate(progress.bands):
                    n_target = band.n_particles + max(band.n_particles // 2, 1)
                    if navigator is None:
                        self._collect_in_band(band, band_idx, n_target)
                    elif not self._continue_size_band(navigator, progress, band_idx, n_target):
                        break
                if sum(band.n_frames for band in progress.bands) == n_frames_before:
                    logger.warning(f"⚠️ No frame left before reaching the precision: {summary}.")
                    return False
        finally:
            self._stats_max_n_frames = None

    @staticmethod
    def _band_frame_descriptors(band: ParticleSizeBand) -> List[Tuple[str, Tuple[float, float], float]]:
        """Frames recorded in a band, as frame-navigator descriptors (label, centre in mm, frame width in mm)."""
        return [
            (frame.frame_id, (frame.center_mm.x, frame.center_mm.y), frame.frame_width_um / 1000.0)
            for frame in band.frames if frame.center_mm is not None
        ]

    def get_particle_stats(self, n_par_target, count_eds_frames: bool = False):
        """
        Analyze frames to collect and save particle size statistics.

        With automatic navigation, the diameter range is split into size bands,
        planned from the smallest particles up (see :func:`plan_size_bands_um`):
        each band is searched with the largest FOV in which its smallest particles
        cover ``powder_meas_cfg.min_par_pixels`` pixels, up to the largest particles
        fitting that FOV. Bands are collected largest first, finer bands in
        subframes of the frames visited.
        The bands, with the number and total area of the frames analysed in each,
        are recorded with the statistics in the sample ledger. The statistics are
        weighted by band, by the particles found per unit area scanned (see
        :func:`autoemx.core.particle_statistics.band_particle_weights`).
        Top-level frames use normal labels (``A0``); finer bands use hierarchical
        ids (``A0_a0``) for sub-frames inside visited parent FOVs.

        The measurement is recorded in the sample ledger
        (``SampleLedger.size_measurements``) with the segmentation of its settings
        (``SampleLedger.segmentations``), and each frame analysed as a frame of the
        ledger with the particles found in it (``SampleLedger.frames``), saved
        regularly. A measurement interrupted with the same target and settings is
        resumed: completed bands are kept, and the interrupted band continues in the
        frames not analysed yet. With different settings, a new measurement starts
        (the interrupted one is kept, with its frames). Frames, measurements and
        statistics of earlier measurements stay in the ledger.

        The statistics (added to ``SampleLedger.size_results`` as the active result)
        use the frames of this measurement and, if
        ``powder_meas_cfg.include_eds_auto_particles``, the frames searched by EDS
        on the sample (or, without them, the particles measured by EDS with
        automatic particle selection).

        With powder_meas_cfg.size_stats_ci_tolerance, the measurement then
        goes on, collecting about 50 % more particles in every band per round,
        until the 95 % confidence intervals of
        powder_meas_cfg.size_stats_ci_quantities are within the tolerance
        (or the frame or time limit is reached, or no frame is left).

        With ``count_eds_frames`` (measurement continuing an EDS acquisition of
        the sample, with automatic navigation), the frames searched by EDS
        (frames of the ledger, if ``powder_meas_cfg.include_eds_auto_particles``)
        count too: the particles
        found in them count towards the target of the size band of their
        diameter, the frames overlapping them are not searched again for the same
        sizes, and finer bands are also searched inside them.

        Parameters
        ----------
        n_par_target : int
            Target number of particles to find **per size band** (particles
            already in the ledger are not counted, unless counted in the frames
            searched by EDS with ``count_eds_frames``); the minimum when
            stopping on precision.
        count_eds_frames : bool, optional
            Whether to count the frames searched by EDS (default: False).

        Returns
        -------
        pandas.DataFrame or None
            Summary statistics DataFrame, or None if no particles were found.
        """
        # Settings of the measurement (size bands change the area limits below)
        powder_meas_cfg = self.powder_meas_cfg.model_copy(deep=True)
        n_par_target = int(n_par_target)
        progress = self._resume_stats_progress(powder_meas_cfg, n_par_target)
        if progress is None:
            progress = self._new_stats_progress(powder_meas_cfg, n_par_target, count_eds_frames)
        self._stats_progress = progress
        self._eds_frames, self._eds_diameters_um = (
            particle_statistics.eds_search_frames(self._stats_eds_scans(progress)) if count_eds_frames
            else ([], np.array([]))
        )

        d_min_um = 2.0 * np.sqrt(self.powder_meas_cfg.min_area_par / np.pi)
        d_max_um = 2.0 * np.sqrt(self.powder_meas_cfg.max_area_par / np.pi)
        min_par_pixels = self.powder_meas_cfg.min_par_pixels

        is_manual = bool(
            getattr(getattr(self.EM, "measurement_cfg", None), "is_manual_navigation", False)
        )
        navigator = getattr(self.EM, "frame_navigator", None)

        if is_manual or navigator is None:
            # Frames chosen by the user: one band, particles of at least min_par_pixels pixels in each frame
            if not progress.bands:
                self._add_stats_band(progress, ParticleSizeBand(
                    d_min_um=d_min_um, d_max_um=d_max_um, n_frames=0, scanned_area_um2=0.0,
                ))
            self._collect_size_band(progress.bands[0], n_par_target)
        else:
            plan = plan_size_bands_um(d_min_um, d_max_um, self._im_width, min_par_pixels, self._frame_width_range_um())
            if not plan or plan[-1][0] > d_min_um * (1 + 1e-9):
                logger.warning(
                    f"⚠️ Particles smaller than {plan[-1][0] if plan else d_max_um:.4g} µm cannot be measured with at "
                    f"least {min_par_pixels} pixels at the smallest frame width of the microscope: they are left out."
                )
            if progress.bands and (
                len(progress.bands) > len(plan)
                or not all(
                    # Same upper bound; the lower one only raised by the microscope (smallest frame it can set)
                    np.isclose(band.d_max_um, planned[1]) and band.d_min_um >= planned[0] * (1 - 1e-9)
                    for band, planned in zip(progress.bands, plan)
                )
            ):
                logger.warning("⚠️ The interrupted measurement used other size bands: a new one starts.")
                if progress.measurement is not None:
                    progress.measurement.status = "interrupted"
                progress = self._new_stats_progress(powder_meas_cfg, n_par_target, count_eds_frames)
                self._stats_progress = progress
            for band_idx, (d_lo, d_hi, fw_planned_um) in enumerate(plan):
                if band_idx < len(progress.bands):
                    if band_idx == len(progress.bands) - 1:
                        # Band interrupted (or last one recorded): continue in the frames not analysed yet
                        band_target = self._band_target(progress.bands[band_idx], band_idx, n_par_target)
                        if not self._continue_size_band(navigator, progress, band_idx, band_target):
                            break
                    continue
                if not self._plan_band_frames(navigator, progress, band_idx, fw_planned_um / 1000.0, reuse_grid=True,
                                              d_range_um=(d_lo, d_hi)):
                    break
                # Frame width set by the microscope (clipped to its range): gives the smallest valid particles
                fw_um = float(
                    navigator.grid_search_fw_mm if navigator.grid_search_fw_mm is not None else fw_planned_um / 1000.0
                ) * 1000.0
                d_pixels_um = min_valid_diameter_um(fw_um / self._im_width, min_par_pixels)
                if d_pixels_um >= d_hi:
                    logger.warning(
                        f"⚠️ Particles smaller than {d_hi:.4g} µm cannot be measured with at least "
                        f"{min_par_pixels} pixels at the frame width set by the microscope: they are left out."
                    )
                    break
                band = ParticleSizeBand(
                    d_min_um=max(d_lo, d_pixels_um), d_max_um=d_hi, frame_width_um=fw_um,
                    n_frames=0, scanned_area_um2=0.0,
                )
                self._add_stats_band(progress, band)
                self._collect_in_band(band, band_idx, self._band_target(band, band_idx, n_par_target))

        precision_reached = None
        if powder_meas_cfg.size_stats_ci_tolerance is not None:
            precision_reached = self._collect_until_precise(progress, navigator if not is_manual else None)

        size_bands = progress.bands
        n_found = sum(band.n_particles for band in size_bands)
        eds_scans = self._stats_eds_scans(progress)
        stats_particles = self._stats_particles(progress, eds_scans)
        # The measurement ends
        self._stats_progress = None
        measurement = progress.measurement
        if measurement is not None:
            measurement.status, measurement.precision_reached = "completed", precision_reached
        if not n_found and not self._eds_diameters_um.size:  # (targets possibly reached by EDS alone)
            logger.error(
                '❌ Could not find any particle. Please check your sample, '
                'or change the contrast/brightness values.'
            )
            self._persist_stats_ledger()
            return None

        result = None
        if self._ledger is not None and measurement is not None:
            result = particle_statistics.recalculate_particle_stats(
                self._ledger, include_eds_auto_particles=powder_meas_cfg.include_eds_auto_particles,
                measurement_id=measurement.measurement_id, segmentation_id=measurement.segmentation_id,
                eds_segmentation_id=self._stats_eds_segmentation(progress),
            )
            if result is not None:
                result = self._ledger.append_size_result(result)
        self._persist_stats_ledger()
        return self.save_particle_statistics(
            particles=stats_particles, bands=size_bands, eds_scans=eds_scans,
            confidence_intervals=result.confidence_intervals_95 if result is not None else None,
        )

    def save_particle_statistics(
        self, output_file_suffix='', particles=None, bands=None, confidence_intervals=None, eds_scans=(),
    ):
        """
        Compute the size statistics of the analysed particles and save the particle sizes, the statistics and the
        histogram in the sample folder (see :func:`autoemx.core.particle_statistics.save_particle_statistics`).

        Parameters
        ----------
        output_file_suffix : str, optional
            String added to output file name
        particles : list of ParticleInfo, optional
            Particles included in the statistics (default: all the analysed particles).
        bands : list of ParticleSizeBand, optional
            Size bands of the measurement, to weight the statistics by band.
        confidence_intervals : dict, optional
            95 % confidence intervals of the statistics (default: computed).
        eds_scans : list of ParticleSizeBand, optional
            Frames searched by EDS on the sample, included in the statistics.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the calculated particle size statistics.
        """
        if os.path.basename(self.results_dir) == cnst.IMAGES_DIR:
            output_dir = os.path.dirname(self.results_dir)
        else:
            output_dir = self.results_dir
        return particle_statistics.save_particle_statistics(
            self.analyzed_pars if particles is None else particles, output_dir, self._sample_id,
            output_file_suffix=output_file_suffix, verbose=self.verbose, bands=bands,
            image_format=self.powder_meas_cfg.size_stats_plot_format,
            confidence_intervals=confidence_intervals,
            eds_scans=eds_scans,
        )

    def _move_and_get_particles_stats_in_frame(self, frame_image=None, pixel_size=None, results_dir=None):
        """
        Move to the next frame and analyze it to detect and characterize particles on the substrate.
    
        This function either acquires a new frame from the electron microscope (if running at the EM)
        or processes a provided image. It detects particles, filters them by size and edge proximity,
        calculates their areas in pixels and micrometers squared, and optionally saves annotated images
        with detected particles and their indices.
    
        Parameters
        ----------
        frame_image : np.ndarray, optional
            Grayscale image of the current frame. If not provided and running at the EM, the image is acquired.
        pixel_size : float, optional
            Pixel size in micrometers. Required if not running at the microscope.
        results_dir : str, optional
            Directory to save result images and masks.
    
        Returns
        -------
        bool or None
            Returns True if at least one valid particle is found in the frame,
            False if no valid particles are detected,
            or None if there are no more frames available (when running at the EM).
    
        Notes
        -----
        - If running at the EM, a new frame is acquired and saved if development_mode=True.
        - If not at the EM, both `frame_image` and `pixel_size` must be provided.
        - The function applies a Gaussian blur to suppress noise before particle detection.
        - Particles are filtered by area and by proximity to the frame edge.
        - The area of each accepted particle is appended to `self.analyzed_pars`.
        - An annotated image with detected particles and their indices is saved for visualization.
        - If no valid particles are found, the mask image is deleted and the function returns False.
        - If at least one particle is found, the microscope focus/contrast/brightness is refreshed if needed.
        """
    
        is_acquired = frame_image is None or pixel_size is None
        if not is_acquired:
            # Image given (e.g. a saved frame reanalysed): used even at the microscope
            self._im_height, self._im_width = frame_image.shape
            self.EM.pixel_size_um = pixel_size
            if results_dir:
                self.results_dir = results_dir
        elif EM_driver.is_microscope_connected():
            self._check_EM_controller_initialization()
            move_to_frame_success = self.EM.go_to_next_frame()
            if move_to_frame_success:
                # Collect image
                frame_image = EM_driver.get_image_data(self._im_width, self._im_height, 1)
                if self.development_mode:
                    cv2.imwrite(os.path.join(self.results_dir, self._sample_id + f'_fr_{self.EM.current_frame_label}.png'), frame_image)
            else:
                # No more frames are available
                return None
        else:
            raise ValueError('This function must be run at the microscope, or it needs to be passed both image and its pixel size')
        
        # Select particles on substrate
        blurred_image = cv2.GaussianBlur(frame_image, (5, 5), 0)
        
        # Get mask of particles on substrate
        save_mask_img = True
        par_mask, mask_img_path = self._get_particles_on_substrate_mask(blurred_image, save_image=save_mask_img)
    
        # Find connected components
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)
        
        # Particles valid for the statistics (within the size limits, not at the frame edge): particle i of the frame
        areas_um2, corrections, marks, counted = self._valid_particles(frame_image, stats, centroids, num_labels)
        par_cntr = len(areas_um2)
        frame_id = self._stats_frame_id() or str(self.EM.current_frame_label)

        # Save the frame images (without annotations, unless annotate_particle_images), and record the frame for the
        # statistics, with the particles to mark on its image. Frames without particles are saved too, so that the
        # measurement can be reanalysed (see autoemx.runners.reanalyze_particle_statistics)
        image_file = self._save_frame_images(
            self._unique_image_stem(f"{self._sample_id}_fr{frame_id}"), frame_image,
            [(x, y, r, i) for i, (x, y, r) in enumerate(marks)],
            mask_image=self.prepare_mask_for_visualization(par_mask),
        )
        self._record_stats_frame(frame_image.shape, frame_id, areas_um2, corrections, marks, image_file,
                                 self._particle_outlines(labels, stats, num_labels, counted))

        # Return True if at least 1 particle was found
        if par_cntr == 0:
            # Mask image written above only in development mode (particle statistics masks are saved as _par_mask)
            if save_mask_img and self.development_mode and os.path.exists(mask_img_path):
                os.remove(mask_img_path)
            return False
        else:
            if is_acquired and not self.EM.measurement_cfg.is_manual_navigation and (time.time() - self.EM.microscope_ctrl._last_EM_adjustment_time > self.EM.microscope_ctrl.refresh_time):
                # Adjust EM focus, contrast and brightness, but only if a particle is actually present
                self.EM.adjust_BCF()
            return True
        