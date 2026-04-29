# DelaunayTriangulation


make clean; make; ./triangle -f inputs/dummy.txt -n 8 -m W -b 1 -i 10

python3 draw.py --points inputs/dummy.txt --triangles outputs/dummy.txt --out outputs/visuals/mesh.png

For divide and conquer, if running the above python3 draw.py command, it is important that you use the sorted points as input:
The algorithm creates a file in inputs called "sorted_{file name}" that you should use as follows:

python3 draw.py --points inputs/sorted_dummy.txt --triangles outputs/dummy.txt --out outputs/visuals/triangulation.png

