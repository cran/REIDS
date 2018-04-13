'''
Created on 6-jul.-2015

@author: lucp8409
'''
#!/usr/bin/env python

if __name__ == '__main__':
	from argparse import ArgumentParser, FileType
	import csv, sys
	arg_parser = ArgumentParser(description='create a CSV file with start '
                                            'line positions and line '
                                            'lengths for a given file')
	arg_parser.add_argument('input_file', type=FileType('r'),
                            help='file to index')
	arg_parser.add_argument('--output_file',
                            help='CSV output file')
	arg_parser.add_argument('--no_header', action='store_true',
                            help='do not add header row to CSV file')
	options = arg_parser.parse_args()
	if options.output_file:
		with open(options.output_file, 'w', newline='') as csvfile:
			csv_writer = csv.writer(csvfile,quoting=csv.QUOTE_NONNUMERIC)
			csv_writer.writerow(['position', 'line_length'])
			c = options.input_file.read(1)
			length = 0
			pos = 0
			line_pos = 0
			while c:
				pos += 1
				if c == '\r':
					pass
				elif c == '\n':
					print(length)
					csv_writer.writerow([line_pos, str(length)])
					length = 0
					line_pos = pos
				else:
					length += 1
				c = options.input_file.read(1)
	else:
		csv_writer = csv.writer(sys.stdout)
    
        
        