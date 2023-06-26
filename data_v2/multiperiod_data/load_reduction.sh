for FILE in *.Pd; do julia load_reduction.jl $FILE 0.7; done
for FILE in *.Qd; do julia load_reduction.jl $FILE 0.7; done