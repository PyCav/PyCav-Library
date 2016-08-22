
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

def green():
    """
    Where would Physicists be without Green functions?

    Parameters
    ----------
    No parameters, just sit back and enjoy the ride.
    """

    import random

    from urllib.request import urlopen
    from bs4 import BeautifulSoup

    url = "https://raw.githubusercontent.com/genneth/dave-green-facts/master/index.html"
    html = urlopen(url)

    soup = BeautifulSoup(html, 'html.parser')
    list = soup.find_all('li')

    print(random.choice(list).text.replace('Dave Green', "\033[0;32mDave Green\033[0m"))
