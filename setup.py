from setuptools import setup, find_packages

setup(
    name="spectro-suite",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python GUI for plotting and comparing 2D spectroscopy datasets.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/spectro-suite", # Replace with your URL
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