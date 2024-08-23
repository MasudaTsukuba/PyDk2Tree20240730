"""
GenerateSampleInputFile20240820.py
2024/8/20, T. Masuda
"""
from datetime import datetime

current_date = datetime.now()
# date_string = current_date.strftime('%Y%m%d%H%M')
date_string = current_date.strftime('%Y%m%d%H')

print(date_string)

number_of_edges = 1000
with open(f'../sample_input_{date_string}.txt', 'w') as outfile:
    for i in range(number_of_edges):
        outfile.write(f'{i+1} {i+2}\n')
