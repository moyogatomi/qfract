from setuptools_rust import RustExtension
from setuptools import setup

setup(
    rust_extensions=[RustExtension("qfract.qfract")],
    name='qfract',
    version="0.0.1",
    description='Rust library with bindings for Python',
    keywords='mandelbrot, coloring, development, experiments, art',
    install_requires=['numpy','pillow','wheel']
)
