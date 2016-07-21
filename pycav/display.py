import matplotlib.pyplot as plt
from IPython.display import HTML
from tempfile import NamedTemporaryFile
import base64
import os


VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""

def anim_to_html(anim,temp,fname,overwrite):
    if temp:
        with NamedTemporaryFile(suffix='.mp4',delete = False) as f:
            anim.save(f.name, fps=20, extra_args=['-vcodec', 'libx264'])
            video = open(f.name, "rb").read()
            os.unlink(f.name)
            assert not os.path.exists(f.name)
    else:
        if not os.path.exists(fname):
            anim.save(fname, fps=20, extra_args=['-vcodec', 'libx264'])
        elif os.path.exists(fname) and overwrite:
        	anim.save(fname, fps=20, extra_args=['-vcodec', 'libx264'])
        video = open(fname, "rb").read()
    encoded_video = base64.b64encode(video).decode('utf-8')
    
    anim._encoded_video = encoded_video

    return anim

def display_animation(fname):
    if hasattr(fname,'_encoded_video'):
    	return HTML(VIDEO_TAG.format(fname._encoded_video))
    else:
        video = open(fname, "rb").read()
        encoded_video = base64.b64encode(video).decode('utf-8')
        return HTML(VIDEO_TAG.format(encoded_video))

def create_animation(anim, temp = False, fname = 'video.mp4', overwrite = True):
    plt.close(anim._fig)
    anim = anim_to_html(anim,temp,fname,overwrite)
    return anim