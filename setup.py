from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    description = fh.read()

setup(
    name="CLJKCov",
    version="0.0.1",
    author="Yuanyuan Zhang",
    author_email="yuanyuan.zhang@noirlab.edu",
    packages=find_packages(), #automatically finds packages in dir
    description="Analyze jackknife covariances from simulations",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/gituser/test-tackage",
    license='MIT',
    python_requires='>=3.8',
    install_requires=['numpy','matplotlib','treecorr','healpy','astropy']
)
