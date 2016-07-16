env = Environment()
env.Replace(CC = "g++")
env.Append(CXXFLAGS = ['-std=c++11'])
Export('env')

common = env.SConscript("BCH_codes/SConscript_object", variant_dir="object/common", duplicate=0)
Export('common')

prog1 = env.SConscript("SConscript_bin", variant_dir="object/BCH_codes", duplicate=0)
test = env.SConscript("bch_tests/SConscript_bin", variant_dir="object/test", duplicate=0)

env.Install('bin', prog1)
env.Install('bin', test)