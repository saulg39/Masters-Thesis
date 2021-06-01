import math
import numpy as np
import os
import xlrd
from xlutils.copy import copy
from xlwt import Workbook

def get_data(row_number):
     THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
     my_file = os.path.join(THIS_FOLDER, 'Data','beam_column_data.xlsx')
     book = xlrd.open_workbook(my_file)
     # get the first worksheet
     sheet = book.sheet_by_index(1)
        # read a row slice


     row = sheet.row_slice(rowx=row_number-1,
                         start_colx=0,
                         end_colx=40)

     Material = row[3].value
     Forming_process = row[5].value
     Bending_type = row[6].value
     Actual_Length = row[7].value
     shape = []
     for i in range(8,15):
          shape.append(row[i].value)
     Material_flat = []
     for i in range(15,21):
          Material_flat.append(row[i].value)
     Material_corner = [row[21].value,[]]
     for i in range(22,28):
          Material_corner[1].append(row[i].value)

     return (Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner)