import matplotlib.pyplot as plt
import numpy as np
import time
import os
import cv2
import warnings
from autoemxsp.tools.utils import EMError


# =============================================================================
# Stage Physical Dimensions
# =============================================================================
'''
Specify stage limits in the microscope's own reference system.
These are also used to transform coordinates from image pixels to absolute microscope stage coordinates.
'''

stage_x_left   = -40  # mm, left limit of stage
stage_x_right  =  40  # mm, right limit of stage
stage_y_top    =  40  # mm, top limit of stage
stage_y_bottom = -40  # mm, bottom limit of stage

'''
Transformation to convert image pixel shifts to stage coordinate shifts.
This encodes axis directionality based on stage limits.
'''
image_to_stage_coords_transform = np.array([
    np.sign(stage_x_right - stage_x_left),   # X axis: +1 if right > left, -1 otherwise
    np.sign(stage_y_bottom - stage_y_top)    # Y axis: +1 if bottom > top, -1 otherwise
])

def frame_rel_to_pixel_coords(pts_rel_coords: np.ndarray, img_width: int, img_height: int) -> np.ndarray:
    """
    Convert relative coordinates (centered at 0, aspect-ratio corrected) to pixel coordinates.
    
    Phenom Coordinate System
    ----------------
    The coordinates are expressed in a normalized, aspect-ratio-correct system centered at the image center:

        - The origin (0, 0) is at the image center.
        - The x-axis is horizontal, increasing to the right, ranging from -0.5 (left) to +0.5 (right).
        - The y-axis is vertical, increasing downward, and scaled by the aspect ratio (height/width):
            * Top edge:    y = -0.5 × (height / width)
            * Bottom edge: y = +0.5 × (height / width)
        
        |        (-0.5, -0.5*height/width)         (0.5, -0.5*height/width)
        |                       +-------------------------+
        |                       |                         |
        |                       |                         |
        |                       |           +(0,0)        |-----> +x
        |                       |                         |
        |                       |                         |
        v  +y                   +-------------------------+
                (-0.5,  0.5*height/width)         (0.5, 0.5*height/width)

    This ensures the coordinate system is always centered and aspect-ratio-correct, regardless of image size.
    
    Parameters
    ----------
    pts_rel_coords : np.ndarray
        Array of shape (N, 2) with (x_rel, y_rel) coordinates.
        - x_rel in [-0.5, 0.5]
        - y_rel in [-0.5 * (H/W), 0.5 * (H/W)]
    img_width : int
        Image width in pixels (W).
    img_height : int
        Image height in pixels (H).

    Returns
    -------
    np.ndarray
        Array of shape (N, 2) with (x_px, y_px) pixel coordinates.
    """
    pts_rel_coords = np.asarray(pts_rel_coords, dtype=float)

    # If it's a single point, reshape to (1, 2)
    if pts_rel_coords.ndim == 1:
        pts_rel_coords = pts_rel_coords[np.newaxis, :]

    aspect_ratio = img_height / img_width

    x_px = (pts_rel_coords[:, 0] + 0.5) * img_width
    y_px = (pts_rel_coords[:, 1] / aspect_ratio + 0.5) * img_height

    return np.column_stack((x_px, y_px))


def frame_pixel_to_rel_coords(pts_pixel_coords: np.ndarray, img_width: int, img_height: int) -> np.ndarray:
    """
    Convert pixel coordinates to relative coordinates (centered at 0, aspect-ratio corrected).
    
    Phenom Coordinate System
    ----------------
    The coordinates are expressed in a normalized, aspect-ratio-correct system centered at the image center:

        - The origin (0, 0) is at the image center.
        - The x-axis is horizontal, increasing to the right, ranging from -0.5 (left) to +0.5 (right).
        - The y-axis is vertical, increasing downward, and scaled by the aspect ratio (height/width):
            * Top edge:    y = -0.5 × (height / width)
            * Bottom edge: y = +0.5 × (height / width)
        
        |        (-0.5, -0.5*height/width)         (0.5, -0.5*height/width)
        |                       +-------------------------+
        |                       |                         |
        |                       |                         |
        |                       |           +(0,0)        |-----> +x
        |                       |                         |
        |                       |                         |
        v  +y                   +-------------------------+
                (-0.5,  0.5*height/width)         (0.5, 0.5*height/width)

    This ensures the coordinate system is always centered and aspect-ratio-correct, regardless of image size.
    
    
    Parameters
    ----------
    pts_pixel_coords : np.ndarray
        Array of shape (N, 2) with (x_px, y_px) pixel coordinates.
    img_width : int
        Image width in pixels (W).
    img_height : int
        Image height in pixels (H).

    Returns
    -------
    np.ndarray
        Array of shape (N, 2) with (x_rel, y_rel) coordinates.
        - x_rel in [-0.5, 0.5]
        - y_rel in [-0.5 * (H/W), 0.5 * (H/W)]
    """
    pts_pixel_coords = np.asarray(pts_pixel_coords, dtype=float)

    # Ensure shape is (N, 2) even if a single point is passed
    if pts_pixel_coords.ndim == 1:
        pts_pixel_coords = pts_pixel_coords[np.newaxis, :]

    aspect_ratio = img_height / img_width  # H/W

    # Map X: [0, W] -> [-0.5, 0.5]
    x_rel = pts_pixel_coords[:, 0] / img_width - 0.5

    # Map Y: [0, H] -> [-0.5*H/W, 0.5*H/W]
    y_rel = (pts_pixel_coords[:, 1] / img_height - 0.5) * aspect_ratio

    return np.column_stack((x_rel, y_rel))

# =============================================================================
# Microscope Navigation Camera Physical Dimensions
# =============================================================================

### Navigation camera
navcam_im_w_mm = 98
"""Width of the navigation camera image, in millimeters."""

stub_w_mm = 12
"""Diameter of the aluminum stubs, in millimeters."""

# Image offsets (if misalignment between navigation camera and SEM is present, otherwise set to 0)
navcam_x_offset = 2
"""X offset between navigation camera and SEM, in pixels."""

navcam_y_offset = 2
"""Y offset between navigation camera and SEM, in pixels."""

# =============================================================================
# SEM acquisition parameters
# =============================================================================
typical_wd = 5 
"""Typical working distance used at SEM, in millimeters."""

im_width = 1920
im_height = 1200
"""SEM image width and height, in pixels."""

# =============================================================================
# Microscope Controls
# =============================================================================
# Connect to electron microscope API
try:
    import PyPhenom as ppi
    phenom = ppi.Phenom()
    acqScanParams = ppi.ScanParams()
    acqScanParams.size = ppi.Size(im_width,im_height)
    acqScanParams.detector = ppi.DetectorMode.All
    acqScanParams.nFrames = 1
    acqScanParams.hdr= False
    acqScanParams.scale = 1.0
    
    is_at_EM = True

except Exception as e:
    is_at_EM = False
    warnings.warn(f'Microscope driver not available: {e}.')


# Limit min magnification to avoid seeing round-shape of lens in image, which can mess with the DQN, as it creates bright points
maxFrameWidth = 0.0011677500002778723 # Corresponds to a magnification of 400x

def standby():
        return phenom.Standby()

def set_electron_detector_mode(detector_name):
    viewingMode = phenom.GetSemViewingMode()
    viewingMode.scanParams.detector = getattr(ppi.DetectorMode, detector_name)
    phenom.SetSemViewingMode(viewingMode)
    
    #TODO add checks that detector names are accepted
    if detector_name != "All":
        raise ValueError(f'Detector type {detector_name} not recognised')

def get_instrument_mode():
    return phenom.GetInstrumentMode()

def get_operational_mode():
    return phenom.GetOperationalMode()

def activate():
    if get_instrument_mode() == ppi.InstrumentMode(2): # Checks if SEM in Standby
        phenom.Activate()
    else:
        pass

def to_SEM(timeout : int = 120):
    start = time.time()
    while True:
        if time.time() - start > timeout:
            raise TimeoutError("SEM did not return to LiveSem within "
                               f"{timeout} s")
        try:
            # Attempt moving to LiveSem OperationalMode
            phenom.MoveToSem()
            break # If the operation succeeds, exit the loop
        except:
            wait_time = 5 #seconds
            print('Phenom busy moving to LiveCamNav.\nWaiting %d seconds before retrying to move to LiveSEM.'% wait_time)
            time.sleep(wait_time) # Wait before trying again
    
    time.sleep(1) # Workaround, otherwise Phenom throws error        



def set_high_tension(voltage):
    phenom.SetSemHighTension(voltage)

def set_beam_current(current):
    phenom.SetSemSpotSize(current)


def get_EDS_analyser_object():
    settings = phenom.LoadPulseProcessorSettings()
    analyzer = ppi.Application.ElementIdentification.EdsJobAnalyzer(phenom)
    analyzer.preset = settings.spot
    
    return analyzer
    


def auto_focus():
    # Pause briefly to ensure autofocus works reliably
    time.sleep(0.2)
    
    # Automatically optimize focus
    phenom.SemAutoFocus()
    
        
    # Get current working distance in mm
    wd = get_current_wd()
    
    return wd

def auto_contrast_brightness():

    # Automatically optimize contrast and brightness
    phenom.SemAutoContrastBrightness()

def adjust_focus(new_wd):
    # Set new working distance
    phenom.SetSemWD(new_wd * 0.001)
    
    
    
def getCurrentPosition():
    pos = phenom.GetStageModeAndPosition().position
    return (pos.x, pos.y)

def moveBy(x, y):

    # Move (relative) position, in meters
    phenom.MoveBy(x*0.001, y*0.001)
    
def move_to(x, y):
    
    # Move (absolute) position, in meters
    phenom.MoveTo(x*0.001, y*0.001)

def get_frame_width():
    current_width = phenom.GetHFW() *1000 #in mm
    return current_width

def convertFrameWtoMag(frameW):
    # Convert to "magnification", for ease of understanding
    mag = ppi.MagnificationFromFieldWidth(frameW)
    # mag = ppi.MagnificationFromFieldWidth(frameW, 0.5) # 0.5 is the displaySize
    return mag

def get_range_frame_width():
    # Probe current mag range, but limit max frame width to minFrameWidth
    range_fw = phenom.GetHFWRange()
    min_fw = range_fw.begin * 1000 # in mm
    max_fw = range_fw.end * 1000 # in mm
    return min_fw, max_fw
    
def set_frame_width(frame_width):
    phenom.SetHFW(frame_width * 0.001)

def zoom(amt):

    current_width = phenom.GetHFW()
    new_width = amt*current_width

    phenom.SetHFW(new_width)

def saveImage(fname='Image.tiff'):
    # Acquire image
    acq = phenom.SemAcquireImage(acqScanParams)
    
    # Check if image already exists, and delete it if it does
    # because ppi.Save() does not overwrite
    if os.path.exists(fname):
        os.remove(fname)
    
    # Save image
    ppi.Save(acq, fname)

def get_image_data(width = 1920, height = 1200, frame_avg = 1):
    #Image must be 8-bit
    acq = phenom.SemAcquireImage(width, height, frame_avg)
    img_data = np.asarray(acq.image)

    return img_data

def showImage():

    img = get_image_data()

    plt.imshow(img, cmap='gray')
    plt.show()

def isFrameWInRange(attempted_frameW):
    
    frameWRange = phenom.GetHFWRange()
    
    if frameWRange.begin <= attempted_frameW <= frameWRange.end:
        validMag = True
    else:
        validMag= False
    
    return validMag

def set_brightness(val):
    phenom.SetSemBrightness(val)
    
def set_contrast(val):
    phenom.SetSemContrast(val)
    
### Set of new commands
def get_current_wd():
    wd = phenom.GetSemWD() # returns value in m
    wd *= 1000 # Convert to mm
    return wd


def acquire_XS_spectral_data(analyzer, x, y, max_collection_time, target_counts, elements = None):
    spectrum_data = None 
    background_data = None
    real_time = None
    live_time = None
    
    # Collect spectrum
    try:
        spotData = analyzer.AddSpot(ppi.Position(x, y), maxTime = max_collection_time, maxCounts = target_counts)
    except Exception as er:
        raise EMError(f"The following error was encountered during X-Ray Spectrum acquisition:\n{er}\n"
        "Spectrum not recorded.")
        return spectrum_data, background_data, real_time, live_time

    # Wait until EDS spectrum collection is over
    analyzer.Wait()
        
    # Get spectrum data
    phenom_spectrum = spotData.spotSpectrum
    
    # Quantify spectrum with Phenom and get background, forcing quantification with elements known to be in the standard
    if elements:
        ppi_elements = [getattr(ppi.Spectroscopy.Element, el) for el in elements]
        try:
            phenom_background = ppi.Spectroscopy.Quantify(phenom_spectrum, ppi_elements).background
        except Exception as er:
            raise EMError(f"The following error was encountered during spectral quantification by EM proprietary software:\n{er}\n"
            "Background not recorded.")
            return spectrum_data, background_data, real_time, live_time
        else:
            background_data = phenom_background.data
    
    spectrum_data = phenom_spectrum.spectrum.data
    real_time = phenom_spectrum.metadata.realTime # Actual time measurement took
    live_time = phenom_spectrum.metadata.liveTime # Detector live time during measurement
    
    return spectrum_data, background_data, real_time, live_time



#### Navigation camera functions


def to_nav():
    # Wake up SEM if necessary
    activate()
     
    # Switch to optical image of sample holder
    try:
        phenom.MoveToNavCam()
    except Exception as e:
        print("Error: Failed to switch to navigation mode:", e)
        return False
    
    return True

def get_navigation_camera_image():
    successful = to_nav()
    if successful:
        # Fix NavCam brightness and contrast for ease of C-tape detection
        phenom.SetNavCamBrightness(0.34)
        phenom.SetNavCamContrast(0.27)
        
        # Set navcam params
        acqCamParams = ppi.CamParams()
        acqCamParams.size = ppi.Size(912, 912)
        acqCamParams.nFrames = 1
    
        # Get navcam image
        acqNavCam = phenom.NavCamAcquireImage(acqCamParams)
        
        # Save and load to convert into numpy format
        temp_f = 'NavCam_temp.tiff'
        ppi.Save(acqNavCam.image, temp_f)
        navcam_im = cv2.imread(temp_f, cv2.IMREAD_COLOR)
        os.remove(temp_f)
        
        return navcam_im
    
    else:
        return None
    
    
    
    