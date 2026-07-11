# Source this in the OOD Desktop terminal to use the shared user-space VirtualGL build.
#   source /projects/bgvl/SummerSchool_2026/DNA/files/legacy/VirtualGL/setup_env.sh
# Then launch GPU-accelerated apps with the EGL back end, e.g.:
#   vglrun -d egl vmd
VGL_HOME=/projects/bgvl/SummerSchool_2026/DNA/files/legacy/VirtualGL
export PATH=$VGL_HOME/bin:$PATH
export LD_LIBRARY_PATH=$VGL_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
