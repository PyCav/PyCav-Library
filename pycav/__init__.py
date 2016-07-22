
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def green():
	import random

	from urllib.request import urlopen
	from bs4 import BeautifulSoup

	url = "http://www.srcf.ucam.org/davegreenfacts/"
	html = urlopen(url)

	soup = BeautifulSoup(html, 'html.parser')
	list = soup.find_all('li')

	print(random.choice(list).text.replace('Dave Green', "\033[0;32mDave Green\033[0m"))
