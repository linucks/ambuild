import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ambuild",
    version="1.0.0",
    author="Jens Thomas",
    author_email="linucks42@gmail.com",
    description="A program for creating polymeric molecular structures.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/linucks/ambuild",
    packages=setuptools.find_packages(),
    install_requires=['numpy'],
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
)
