import nibabel as nib
import matplotlib.pyplot as plt

import os



for p in range(10,90,1):
    
    print(p)

    dirname = os.path.dirname(__file__)
    path = os.path.join(dirname, 'training', 'patient0{0}'.format(p), 'patient0{0}_frame01_gt.nii.gz'.format(p))
    #path = os.path.join(dirname, 'training', 'patient00{0}'.format(p), 'patient00{0}_frame01.nii.gz'.format(p))

    # Change the path to your path
    my_img  = nib.load(path)
    nii_data = my_img.get_fdata()
    nii_aff  = my_img.affine
    nii_hdr  = my_img.header
    print(nii_aff ,'\n',nii_hdr)
    print("nii_data.shape -- ")
    print(nii_data.shape)
    plt_i = 0

    if(len(nii_data.shape)==3):
        for slice_Number in range(nii_data.shape[2]):
            plt.axis('off')
            plt.imshow(nii_data[:,:,slice_Number ])
            plt.savefig(dirname + 'patient{0}-gt-fig-{1}.png'.format(p, plt_i+1))
            #plt.savefig(dirname + 'patient{0}-pred-fig-{1}.png'.format(p, plt_i+1), transparent=True)
            #plt.show()
            plt_i = plt_i + 1
        
    if(len(nii_data.shape)==4):
        for frame in range(nii_data.shape[3]):
            for slice_Number in range(nii_data.shape[2]):
                plt.axis('off')
                plt.imshow(nii_data[:,:,slice_Number,frame])
                plt.savefig(dirname + 'patient{0}-gt-fig-{1}.png'.format(p, plt_i+1))
                #plt.savefig(dirname + 'patient{0}-pred-fig-{1}.png'.format(p, plt_i+1), transparent=True)
                #plt.show()
                plt_i = plt_i + 1
            
            
            

    