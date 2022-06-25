#!/usr/local/env python3

import os
import gzip
from concurrent import futures
from argparse import ArgumentParser
import multiprocessing as mp
from itertools import repeat, islice


__author__ = 'duceppemo'
__version__ = '0.1'


class FastqParser(object):
    @staticmethod
    def make_chunks(file_handle, size):
        while True:
            chunk = list(islice(file_handle, size))
            if not chunk:
                break
            yield chunk

    @staticmethod
    def read_and_filter_entry(fastq_chunk, range_start, range_stop):
        """
        Reads a single fastq entry and filter it based on 'ch' value
        :param fastq_chunk: list. Length of 4 [header, seq, extra, qual].
        :param range_start: int.
        :param range_stop: int.
        :return: dict. { sample_name: [header, seq, extra, qual] }
        """
        lines = list()
        fastq_dict = dict()
        for line in fastq_chunk:
            line = line.rstrip()
            lines.append(line)
            if len(lines) == 4:
                # header, seq, extra, qual = single_fastq_entry
                # @b13de425-8e68-4cb7-a946-e524ac5ae492 runid=97039b6eb894fd91e4dcbd537db29268f9ca15a7 sampleid=tets1_R10 read=45559 ch=270 start_time=2022-06-16T00:20:20Z model_version_id=2021-11-17_dna_r10.4_minion_promethion_1024_67af0493 barcode=barcode07
                seq_id = lines[0].split()[0]
                channel = [x for x in lines[0].split() if 'ch=' in x][0].split('=')[1]
                channel = int(channel)
                if channel in range(range_start, range_stop + 1):
                    fastq_dict[seq_id] = lines
        return fastq_dict

    @staticmethod
    def iterate_fastq_parallel(input_fastq, output_folder, range_start, range_stop, cpu, parallel):
        """
        Split the fastq files into chunks of 1000 entries (4000 lines) for parallel processing.
        :param input_fastq:
        :param output_folder:
        :param range_start:
        :param range_stop:
        :param cpu:
        :param parallel:
        :return:
        """
        # Name
        sample_name = os.path.basename(input_fastq).split('.')[0].split('_')[0]
        sample_name = sample_name.replace('_pass', '')
        sample_name = sample_name.replace('_fail', '')

        # Chunk fastq files and run chunks in parallel
        fastq_dict = dict()
        with gzip.open(input_fastq, "rt") if input_fastq.endswith('.gz') else open(input_fastq, "r") as f:
            pool = mp.Pool(int(cpu / parallel))
            jobs = [pool.apply_async(FastqParser.read_and_filter_entry, [chunk, range_start, range_stop])
                    for chunk in FastqParser.make_chunks(f, 4000)]
            results = [j.get() for j in jobs]
            pool.close()
            pool.join()
            pool.terminate()  # Needed to do proper garbage collection?

            # Update self.sample_dict with results from every chunk
            for d in results:
                fastq_dict.update(d)  # Do the merge

        # Write filtered fastq to file
        output_fastq = '{}/{}.fastq.gz'.format(output_folder, sample_name)
        with gzip.open(output_fastq, 'wb') as f:
            for seq_id, fastq_entry_list in fastq_dict.items():
                f.write('{}\n'.format('\n'.join(fastq_entry_list)).encode('ascii'))

    @staticmethod
    def parallel_process_fastq(fastq_list, output_folder, range_start, range_stop, cpu, parallel):
        """
        Parallel process the fastq files in the input folder
        :param fastq_list:
        :param output_folder:
        :param range_start:
        :param range_stop:
        :param cpu:
        :param parallel:
        :return:
        """
        with futures.ProcessPoolExecutor(max_workers=parallel) as executor:
            for results in executor.map(FastqParser.iterate_fastq_parallel, fastq_list, repeat(output_folder),
                                        repeat(range_start), repeat(range_stop), repeat(cpu), repeat(parallel)):
                pass
            # args = ((fastq, output_folder, range_start, range_stop, int(cpu/parallel), parallel)
            #         for fastq in fastq_list)
            # for results in executor.map(lambda p: FastqParser.iterate_fastq_parallel(*p), args):
            #     pass


class PoreRangeSplitter(object):
    def __init__(self, args):
        """Define objects based on supplied arguments"""
        self.args = args
        self.input_folder = args.input
        self.output_folder = args.output
        self.pore_range = args.range
        self.cpu = args.threads
        self.parallel = args.parallel

        # run the script
        self.run()

    def run(self):
        # Check range
        range_start, range_stop = PoreRangeSplitter. check_range(self.pore_range)

        # Process fastq
        fastq_list = PoreRangeSplitter.list_fastq_files(self.input_folder)
        FastqParser.parallel_process_fastq(fastq_list, self.output_folder, range_start, range_stop,
                                           self.cpu, self.parallel)

    @staticmethod
    def check_range(my_range):
        try:
            range_start, range_stop = my_range.split('-')
        except ValueError:
            raise Exception('Range format needs to be "start-end".')

        # Convert to integer
        range_start = int(range_start)
        range_stop = int(range_stop)

        # Check values
        if range_start > range_stop:
            raise Exception('The first value of the range must be lower than the second one.')

        # if (range_start or range_stop < 0) or (range_start or range_stop > 512):
        if 0 > (range_start or range_stop) > 512:
            raise Exception('Pore range values must be between 0-512.')

        return range_start, range_stop

    @staticmethod
    def list_fastq_files(input_folder):
        fastq_list = list()

        # Walk input folder recursively to find fastq files
        if input_folder and os.path.isdir(input_folder):
            for root, directories, filenames in os.walk(input_folder):
                for filename in filenames:
                    absolute_path = os.path.join(root, filename)
                    if os.path.isfile(absolute_path) and filename.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                        fastq_list.append(absolute_path)

        # Check if input_fastq_list is not empty
        if not fastq_list:
            raise Exception("No fastq file found in %s!" % input_folder)

        return fastq_list


if __name__ == '__main__':
    max_cpu = mp.cpu_count()

    parser = ArgumentParser(description='Split Nanopore reads based on channel number.')
    parser.add_argument('-i', '--input', metavar='/path/to/folder/with/fastq',
                        required=True, type=str,
                        help='Input folder with fastq file(s), gzipped or not.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder',
                        required=True, type=str,
                        help='Output folder.')
    parser.add_argument('-t', '--threads', metavar='{}'.format(max_cpu),
                        required=False, type=int, default=max_cpu,
                        help='Number of CPU. Default {}'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='{}'.format(max_cpu),
                        required=False, type=int, default=max_cpu,
                        help='Number of samples to process in parallel. Default 4.')
    parser.add_argument('-r', '--range', metavar='0-256',
                        required=True, type=str,
                        help='Pore range to keep. Range between 0-512')

    # Get the arguments into an object
    arguments = parser.parse_args()

    PoreRangeSplitter(arguments)
