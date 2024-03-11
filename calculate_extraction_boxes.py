# Alter these values:
pixel_sizes = [1.14, 1.08]
minimum_diameter = 490 # in Angstrom

# than run the rest of the code

# good box sizes (you might alter these as well) (retrieved from: https://blake.bcm.edu/emanwiki/EMAN2/BoxSize)
good_boxsizes = [24, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 84, 96, 100, 104, 112, 120, 128, 132, 140, 168, 180, 192, 196, 208, 216, 220, 224, 240, 256, 260, 288, 300, 320, 352, 360, 384, 416, 432, 440, 448, 480, 512, 540, 560, 576, 588, 600, 630, 640, 648, 672, 686, 700, 720, 750, 756, 768, 784, 800, 810, 840, 864, 882, 896, 900, 960, 972, 980, 1000, 1008, 1024 ]

import numpy as np

# some functions
is_even = lambda num: (num % 2) == 0
is_within_threshold = lambda diff: abs(diff) <= 1e-4
closest = lambda target: min(good_boxsizes, key=lambda x: abs(x - target))
previous_box = lambda value: good_boxsizes[good_boxsizes.index(value) - 1]

def get_downsampled_ps(initial_pixel_size, target_box=downsample_to, target_df=1, return_extraction=False):
    target_dpx = target_df / initial_pixel_size
    ext_initial = np.floor(target_dpx * target_box)

    if not is_even(ext_initial):
        ext_initial += 1

    downsampled_factor = ext_initial/target_box
    downsampled_ps = initial_pixel_size*downsampled_factor 
    
    if return_extraction:
        return ext_initial
        
    return downsampled_ps

def ps_differences(lst):

    arr = np.array(lst)
    differences = arr[:, np.newaxis] - arr
    differences = differences[np.triu_indices(len(arr), k=1)]

    return differences

target_dfs = np.arange(max(pixel_sizes),max(pixel_sizes)+0.5, 0.0005)
success = False
target_box = closest(minimum_diameter / max(pixel_sizes) )
debug_counter = 0
while success == False:
    for n in target_dfs:
        downsampled_ps_collected = [get_downsampled_ps(initial_pixel_size=k, target_box=target_box, target_df=n) for k in pixel_sizes]
        diff = ps_differences(downsampled_ps_collected)

        criteria = [is_within_threshold(k) for k in diff]

        if all(criteria):
            extraction_collected = {k : get_downsampled_ps(initial_pixel_size=k, target_df=n, target_box=target_box, return_extraction=True) for k in pixel_sizes}
            success = True
            break
        else:
            continue
    
    debug_counter += 1
    
    # try next smaller boxsize for downsampling
    if not success:

        new_target_box = previous_box(target_box)

        if new_target_box > target_box:
            break
        else:
            target_box = new_target_box
    
if success:
    target_pixelsize = "{:.5f}".format(n)
    print(f'Target pixelsize: {target_pixelsize}')
    print(f'Downsample the extracion to {target_box} pixel')
    for px, box in extraction_collected.items():
        angstrom = "{:.2f} Å".format(px * box)
        print(f'Extract the particles at pixelsize of {px} with a box of {int(box)} pixel ({angstrom})')
    
    formatted_differences = "\n\t".join([str(d) for d in diff])
    
    print(f'Pixelsize variations: \n\t{formatted_differences}')
else:
    print('No combination found.')
