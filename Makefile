all: test python

test:
	g++ -std=c++2a -I. -I/usr/include test_sph.cpp

python:
	g++ sph_py.cpp -I. -I/usr/include `python -m pybind11 --includes` -shared -fPIC -o cppsph`python3-config --extension-suffix`
