import matplotlib.pyplot as plt
from IPython.display import HTML
import base64
import os

VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""

def anim_to_html(anim,temp,fname,overwrite):
    if temp:
        anim._html5_video = anim.to_html5_video()
    else:
        _fps = int(anim.save_count/anim._interval)
        if _fps == 0:
            _fps = 20
        if not os.path.exists(fname):
            anim.save(fname, fps=_fps, extra_args=['-vcodec', 'libx264'])
        elif os.path.exists(fname) and overwrite:
            anim.save(fname, fps=_fps, extra_args=['-vcodec', 'libx264'])
        video = open(fname, "rb").read()

        encoded_video = base64.b64encode(video).decode('utf-8')
        anim._encoded_video = encoded_video

    return anim

def display_animation(fname):
    if hasattr(fname,'_encoded_video'):
    	return HTML(VIDEO_TAG.format(fname._encoded_video))
    elif hasattr(fname,'_html5_video'):
        return HTML(fname._html5_video)
    else:
        video = open(fname, "rb").read()
        encoded_video = base64.b64encode(video).decode('utf-8')
        return HTML(VIDEO_TAG.format(encoded_video))

def create_animation(anim, temp = False, fname = 'video.mp4', overwrite = True):
    plt.close(anim._fig)
    anim = anim_to_html(anim,temp,fname,overwrite)
    return anim
