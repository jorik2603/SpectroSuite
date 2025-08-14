from setuptools import setup, find_packages

setup(
    name="SpectroSuite",
    version="0.1.0",
    author="Jorik Schaap",
    author_email="j.h.schaap@utwente.nl",
    description="A Python GUI for plotting and comparing 2D spectroscopy datasets.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jorik2603/SpectroSuite", # Replace with your URL
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "cmcrameri",
    ],
    python_requires='>=3.8',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    entry_points={
        "console_scripts": [
            "spectrosuite=spectro_suite.__main__:main",
        ],
    },
)