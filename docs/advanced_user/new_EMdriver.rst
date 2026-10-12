.. _advanced_new_EM_driver:

===================================
Add new Electron Microscope Driver
===================================

This guide explains how to write the driver that lets `AutoEMX` control a new microscope through its API.
A driver is only needed for automated acquisition: fitting, quantification and analysis of spectra work without it.

Before you begin, check whether a driver for your microscope is already available in
`autoemx/microscope_drivers <https://github.com/CederGroupHub/AutoEMX/tree/main/autoemx/microscope_drivers>`_.

Start from the Phenom XL driver
-------------------------------

The driver of the Phenom XL, ``autoemx/microscope_drivers/PhenomXL.py``, is the template: each function is
documented with its parameters, units and returned values, and the microscope is controlled only through calls to
its API, ``PyPhenom``.

1. Copy ``PhenomXL.py`` to ``autoemx/microscope_drivers/<microscope_ID>.py``.
2. Replace the ``PyPhenom`` calls with the equivalent calls of your microscope API, keeping the names, signatures,
   units and returned values of the functions. Keep raising ``EMError`` when an operation fails.
3. Set the constants at the top of the file for your microscope: stage limits, navigation camera, image size,
   typical working distance, and ``RAW_SPECTRUM_EXT``, the format of the spectra exported by your EDS software.
4. Add the calibrations of the microscope in ``autoemx/calibrations/<microscope_ID>/`` (see
   :ref:`advanced_new_XSp_calibs`).

The driver is found from the microscope ID set in the configurations (``MicroscopeConfig.ID``): `AutoEMX` loads
``microscope_drivers/<microscope_ID>.py``. Instruments of the same model, each with its own calibration folder, can
share one driver: set ``"driver"`` in their ``detector_model.json`` (e.g. ``"driver": "PhenomXL"``).

What `AutoEMX` uses from the driver
-----------------------------------

`AutoEMX` calls the driver functions by name, so all of the following must be defined.

- **Connection**: ``connect_to_microscope``, ``is_microscope_connected``, and the flag ``is_at_EM``. Import the API
  inside ``connect_to_microscope``, so that the driver can be loaded on computers without it.
- **Microscope state**: ``activate``, ``standby``, ``to_SEM``, ``to_nav``.
- **Beam and detectors**: ``set_high_tension``, ``set_beam_current``, ``set_electron_detector_mode``.
- **Stage**: ``move_to``, and the stage limits ``stage_x_left``, ``stage_x_right``, ``stage_y_top``,
  ``stage_y_bottom`` with ``image_to_stage_coords_transform`` (axis directions).
- **Imaging**: ``get_image_data``, ``get_frame_width``, ``set_frame_width``, ``get_range_frame_width``,
  ``auto_focus``, ``adjust_focus``, ``get_current_wd``, ``auto_contrast_brightness``, ``set_brightness``,
  ``set_contrast``, and the constants ``im_width``, ``im_height``, ``typical_wd``.
- **Image coordinates**: ``frame_rel_to_pixel_coords`` and ``frame_pixel_to_rel_coords``, converting between pixels
  and the coordinates used by the API (see the description of the Phenom coordinate system in the driver).
- **Navigation camera**: ``get_navigation_camera_image``, and its calibration ``navcam_im_w_mm``,
  ``navcam_x_offset``, ``navcam_y_offset``.
- **EDS**: ``get_EDS_analyser_object`` and ``acquire_XS_spectral_data``, which returns the spectrum, the background
  quantified by the EDS software (or ``None``), the real time and the live time.

Optional:

- ``get_position``: current stage position (mm), read in manual navigation after you move the stage. Without it,
  `AutoEMX` keeps using the last known position (initially the sample center).

Testing the driver
------------------

``tests/Test_EM_driver.py`` runs each driver function once with safe example values, for manual verification on
the microscope: set ``microscope_ID``, uncomment the functions to test, and check the result of each one.

Consider contributing your driver to `AutoEMX`, so that it becomes available to other users of the same microscope.
