#env = Environment()
Import('env')


cpp = ["/usr/local/Cellar/glog/0.3.4/lib", \
       "/Users/tarunchitra/Dropbox/Dimers/blosc", \
       "/anaconda/pkgs/hdf5-1.8.17-1/lib" ]

flg = ["-I/anaconda/pkg/hdf5-1.8.17-1/include",\
       "-I/usr/local/Cellar/glog/0.3.4/include", \
       "-I/usr/local/Cellar/gflags/2.1.2/include", \
       "-I/anaconda/pkgs/hdf5-1.8.17-1/include", \
       "-I/Users/tarunchitra/Dropbox/Dimers/blosc" ]

env.Replace(CPPPATH=cpp)
env.Append(CFLAGS=flg, CXXFLAGS=flg)

cflags = ['-Wall', '-Werror', '-std=c++11', '-O0'] # '-g', '-ggdb', '-O0' ]
# ndebug
cflags.append("-DNDEBUG")
# Debug: 
#cflags.append("-D__DEBUG")
lflags = ['-lglog', '-lhdf5' ]
#lflags = ['-lhdf5']

env.Append(
        LINKFLAGS = lflags,
        CXXFLAGS = cflags,
        )

env.SConscript('src/SConscript', 'env')
