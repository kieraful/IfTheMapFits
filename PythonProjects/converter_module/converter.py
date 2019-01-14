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


def csv_to_pcd(files, outfilename='All_points.pcd'):

    pcl_file_tmp = open('All_points_temp.pcd', 'a')
    num_points = 0
    count = 0

    for file in files:
        count += 1
        print('\tWorking on file:\t{}\t{}\r'.format(count, num_points))
        csv_file = open(file)

        csv_reader = csv.reader(csv_file, delimiter=',')

        linecount = 0
        for line in csv_reader:
            linecount += 1
            if linecount != 1:
                num_points += 1
                pcl_file_tmp.write('{X}\t{Y}\t{Z}\t{INTENSITY}\n'.format(X=line[3],
                                                                     Y=line[4],
                                                                     Z=line[5],
                                                                     INTENSITY=line[6]))

    pcl_file_tmp.close()
    pcl_file = open(outfilename, 'a')
    # Write the header for PCD files
    write_pcl_header(pcl_file, num_points)
    # Write the points
    with open('All_points_temp.pcd') as infile:
        for line in infile:
            pcl_file.write(line)
    pcl_file.close()
    # Remove temp file
    os.remove('All_points_temp.pcd')

def convert_dir(path, data_type='csv', output='txt'):
    os.chdir(path)
    count = 0
    print('\nLooking for {} files...\n\n'.format(data_type))

    if output == 'pcd':
        csv_to_pcd(glob.glob("*.{}".format(data_type)))
    else:
        for data_file in glob.glob("*.{}".format(data_type)):
            count += 1
            print('\tWorking on file:\t{}\r'.format(count))
            if output == 'txt':
                csv_to_txt(data_file)
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

