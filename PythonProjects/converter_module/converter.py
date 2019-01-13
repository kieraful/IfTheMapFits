"""Python Converter Utility for If the Map Fits
"""


# Python Default Packages
# import math
import csv
import os
import glob



def csv_to_txt(filename, outfilename='All_points.txt'):

    csv_file = open(filename)
    text_file = open(outfilename, 'a')

    csv_reader = csv.reader(csv_file, delimiter=',')

    linecount = 0
    for line in csv_reader:
        linecount += 1
        if linecount != 1:
            text_file.write('{TIME}\n{X}\t{Y}\t{Z}\t{INTENSITY}\n'.format(TIME=line[10],
                                                                          X=line[3],
                                                                          Y=line[4],
                                                                          Z=line[5],
                                                                          INTENSITY=line[6]))

def csv_to_pcd(filename, outfilename='All_points.pcd'):

    num_points = find_num_points(filename);
    csv_file = open(filename)
    pcl_file = open(outfilename, 'a')

    csv_reader = csv.reader(csv_file, delimiter=',')

    # Write the header for PCD files
    write_pcl_header(pcl_file, num_points )

    linecount = 0
    for line in csv_reader:
        linecount += 1
        if linecount != 1:
            pcl_file.write('{X}\t{Y}\t{Z}\t{INTENSITY}\n'.format(X=line[3],
                                                                 Y=line[4],
                                                                 Z=line[5],
                                                                 INTENSITY=line[6]))
def convert_dir(path, data_type='csv', output='txt'):
    os.chdir(path)
    count = 0
    print('\nLooking for {} files...\n\n'.format(data_type))
    for data_file in glob.glob("*.{}".format(data_type)):
        count += 1
        print('\tWorking on file:\t{}\r'.format(count))
        if output == 'txt':
            csv_to_txt(data_file)
        elif output == 'pcd':
            csv_to_pcd(data_file)
    print('\n\nProcessing finished. {} files converted to Point Cloud {} format'.format(count, output))

def write_pcl_header(file, n_points):

    file.write('VERSION .7\n')
    file.write('FIELDS x y z int\n')
    file.write('SIZE 4 4 4 4\n')
    file.write('TYPE F F F F\n')
    file.write('COUNT 1 1 1 1\n')
    file.write('WIDTH {}\n'.format(n_points))
    file.write('HEIGHT 1\n')
    file.write('VIEWPOINT 0 0 0 1 0 0 0\n')
    file.write('POINTS {}\n'.format(n_points))
    file.write('DATA ascii\n')


def find_num_points(csv_filename):
    csv_file = open(csv_filename)
    num_points = sum(1 for line in csv.reader(csv_file, delimiter=','))
    csv_file.close()
    return num_points

