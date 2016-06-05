#!/bin/bash

head -1 matrix.reorder > matrix.reorder.ccs | tail -n +2 matrix.reorder | sort -n -k 2,2 -k 1,1 >> matrix.reorder.ccs
head -1 matrix.reorder > matrix.reorder.crs | tail -n +2 matrix.reorder | sort -n -k 1,1 -k 2,2 >> matrix.reorder.crs


