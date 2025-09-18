# Create animation from a series of images
# LJK
# 

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as mpimg
import os

def create_animation(image_folder, output_path, interval=200):
    fig, ax = plt.subplots(1,1,figsize=(14,10))
    ax.set_axis_off()
    
    images = []
    img_files = sorted([f for f in os.listdir(image_folder) if f.endswith('.png')])

    for img_file in img_files:
        img_path = os.path.join(image_folder, img_file)
        img = mpimg.imread(img_path)
        img_plot = ax.imshow(img, animated=True)
        images.append([img_plot])
            
    ani = animation.ArtistAnimation(fig, images, interval=interval, blit=True, repeat_delay=1000)
    ani.save(output_path,dpi=350)
    plt.close()

image_folder = data_dir + "G4_traj_paper_figs/traj_animation/"
output_path = "traj_anim.mp4"
create_animation(image_folder, image_folder+output_path, 400)
