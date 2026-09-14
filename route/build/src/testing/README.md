# Unit testing 
This directory contains tests of individual, isolated software components, such as subroutines, functions, or classes, to verify they work exactly as expected.

Each test is stored in a subdirectory named `<test_name>` containing a cmake file (`CMakeLists.txt`) and `<test_name>/src/*.f90`. 

The components to be tested should be called directly from the software source code located under: 

* `route/build/src`
* `route/build/src/standalone`
* `route/build/cpl`

To build an executable for each test:

```bash
cd <test_name>
cmake -S . -B cmake_build
cmake --build cmake_build
```

To run the test,
```bash
./cmake_build/test.exe
```
