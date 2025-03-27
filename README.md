# AlignmentActive

The underlying [AlignmentRepaC repository](https://github.com/caiks/AlignmentRepaC) is a fast C++ implementation of some of the *practicable inducers* described in the paper *The Theory and Practice of Induction by Alignment* at https://greenlake.co.uk/. The AlignmentRepaC repository depends on the [AlignmentC repository](https://github.com/caiks/AlignmentC) for the underlying *model* framework. 

The AlignmentActive repository brings together a *history* and a *model* to define a thread-safe structure for realtime *aligned induction*. The AlignmentActive repository depends on the AlignmentRepaC repository.

## Download and build

The `AlignmentActive` module requires [modern C++](https://en.cppreference.com/w/) version 17 or later to be installed.

For example, in Ubuntu bionic (18.04),
```
sudo apt-get update -y && sudo apt install -y git g++ cmake

```
Then download the zip file or use git to get the underlying rapidjson, AlignmentC and AlignmentRepaC repositories, and the AlignmentActive repository -
```
git clone https://github.com/Tencent/rapidjson.git
git clone https://github.com/caiks/AlignmentC.git
git clone https://github.com/caiks/AlignmentRepaC.git
git clone https://github.com/caiks/AlignmentActive.git

```
Note that the AlignmentActive repository has more than one version. To clone another branch, e.g.
```
git clone -b v01 --single-branch https://github.com/caiks/AlignmentActive.git

```
Then build -
```
cd
mkdir -p AlignmentActive_build
cd ~/AlignmentActive_build
cmake -DCMAKE_BUILD_TYPE=RELEASE ../AlignmentActive
make

```


