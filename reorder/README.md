cd src
make
./ccs
./ccs_reorder
./bcrs

The Mflops will be reported to the standard output, followed with verification message.

matrix.reorder stores the reordered matrix in AIJ format.
matrix.reorder.ccs stores the reorder matrix in AIJ format and sorted based on the column for ccs usage.
matrix.reorder.crs stores the reorder matrix in AIJ format and sorted based on the row for crs usage.

make clean

make pristine

