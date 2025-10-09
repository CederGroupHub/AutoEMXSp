#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:34:01 2025

@author: Andrea
"""

#%% Test particle finding in frame

# frame_images_path = '/Users/Andrea_1/Desktop/Work/SEM EDX automation/Manual tests EDX/SEM images for particle finder/Frame images/'
# frame_images_path = '/Users/Andrea_1/Desktop/Work/SEM EDX automation/Manual tests EDX/David particle counting'

# sample_ID = '15_8'

# frame_path = os.path.join(frame_images_path, sample_ID + '.tiff')

# # with tifffile.TiffFile(frame_path) as tif:
# #     # Print TIFF tags (metadata)
# #     for page in tif.pages:
# #         for tag in page.tags.values():
# #             tag_name = tag.name
# #             tag_value = tag.value
# #             print(f'Tag: {tag_name}, Value: {tag_value}')

# pixel_size_frame_um = 254 / 1920
# pixel_size_frame_um

# frame_image = cv2.imread(frame_path, cv2.IMREAD_GRAYSCALE)
# frame_image.shape
# # cv2.imshow('image',frame_image)
# # Remove databar
# frame_image[1200:,:] = 0

# # # Display the image in a window
# # cv2.imshow('', frame_image)

# particle_finder = SEM_Particle_Finder('Test', (0,0), par_brightness_thresh = 100,
#                                       max_area_par = max_area_par, min_area_par = min_area_par, verbose = True)
# results_dir = frame_images_path
# particle_finder.par_areas = []
# particle_finder.sample_ID = sample_ID
# particle_finder.frame_cntr = 0
# particle_finder.current_pos = np.array([0,0])
# particle_finder._get_particles_stats_in_frame(frame_image, pixel_size_frame_um,  results_dir)