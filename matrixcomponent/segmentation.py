#!/usr/bin/env python3
"""
Tasks
Get example output for tests - Ted
Parser for ODGI Bin file format - Ted
Component Segmentation Detection - Josiah and Joerg
  Python memory object model - Josiah
Output format
"""
from collections import defaultdict
from typing import List, Tuple, Set, Dict

from nested_dict import nested_dict

from matrixcomponent.matrix import Path, PangenomeSchematic, Component, LinkColumn
import os
import logging
import argparse
import matrixcomponent

import matrixcomponent.JSONparser as JSONparser

LOGGER = logging.getLogger(__name__)
"""logging.Logger: The logger for this module"""


def segment_matrix(matrix: List[Path]) -> PangenomeSchematic:
    print(f"Starting Segmentation process on {len(matrix)} Paths.")
    schematic = PangenomeSchematic([])
    incoming, outgoing, dividers = find_dividers(matrix)
    start_pos = 0
    for valid_start in sorted(list(dividers)):
        if valid_start != 0:
            schematic.components.append(Component(start_pos, valid_start - 1))
        start_pos = valid_start

    # TODO: populate all link columns onto schematic
    for component in schematic.components:
        # TODO: order columns based on traversal patterns,
        # TODO: insert additional columns for higher copy number
        for arriving_pos, participants in outgoing[component.last_bin].items():
            leaving = LinkColumn(component.last_bin,
                                 arriving_pos,
                                 participants=participants)
            component.departures.append(leaving)
        for origin_pos, participants in incoming[component.first_bin].items():
            entering = LinkColumn(origin_pos,
                                  component.first_bin,
                                  participants=participants)
            component.arrivals.append(entering)

    return schematic


def find_dividers(matrix: List[Path]) -> Tuple[Dict[int, Dict[int, set]],
                                               Dict[int, Dict[int, set]], Set[int]]:
    max_bin = 13176+1 # TODO: find this
    leaving = nested_dict(2, set)  # containing the set of participating Paths on that link column
    entering = nested_dict(2, set)  # list of indices of new components
    dividers = {0}  # all start positions of components
    for path in matrix:
        print(f"Segmenting {len(path.bins)}")
        for link in path.links:  # Links are generated by odgi based
            upstream, downstream = link.upstream, link.downstream
            # Is the gap range anywhere else in this individual?
            # What if downstream < upstream?
            divider_verified = downstream < upstream
            if not divider_verified:
                missing_range = list(range(upstream + 1, downstream))
                for i in missing_range:
                    if i in path:
                        divider_verified = True
                        break  # stop as soon as we have confirmation
            if divider_verified:
                # if (upstream + 1) in leaving.keys() :
                #     print(f"Found inherited rearrangement {upstream+1}")
                leaving[upstream][downstream].add(path)  # the first position of the new component
                entering[downstream][upstream].add(path)
                dividers.add(upstream + 1)
                dividers.add(downstream)
                # TODO: insert prevarications about exact position
                # Divider should be somewhere in here
                # Tolerable range?
                # Stack up others using the same LinkColumn
    return entering, leaving, dividers


def discard_useless_links(matrix: List[Path]):
    """https://github.com/vgteam/odgi/issues/48
    Links that simply span a gap in the matrix can be discarded"""
    for path in matrix:
        keep = []
        for link in path.links:  # Links are generated by odgi based
            missing_range = list(range(link.upstream + 1, link.downstream))
            # Is the gap range anywhere else in this individual?
            if any([i in path for i in missing_range if i > 0]):
                keep.append(link)
        path.links = keep  # all other Paths get deleted


def setup_logging(output_dir):
    """Setup the logging, add a log file"""
    log_name = os.path.join(output_dir, 'log')
    handler = logging.FileHandler(log_name)
    handler.setLevel(args.log_level)
    handler.setFormatter(logging.Formatter(matrixcomponent.LOGGING_FORMAT_STR,
                                           datefmt=matrixcomponent.LOGGING_DATE_FORMAT))
    logging.getLogger().addHandler(handler)


# Helper class to allow multi-line help messages for argparse user parameters:
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_arguments():
    """Create the command line interface and return the command line arguments

    Returns
    -------
    Namespace
        The command line arguments

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="")

    parser.add_argument('-j', '--json-file',
                            dest='json_file',
                            required=True,
                            help='input JSON file')

    parser.add_argument('-o', '--out-folder',
                        dest='output_folder',
                        required=True,
                        help='output folder')

    parser.add_argument('-l', '--log-level',
                        default='DEBUG',
                        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
                        help='level of logging verbosity. DEBUG is most verbose')

    args = parser.parse_args()

    return args

def main():
    global args
    args = get_arguments()
    setup_logging(args.output_folder)
    LOGGER.info("starting...\n")
    Paths = JSONparser.parse(args.json_file)
    schematic = segment_matrix(Paths)


if __name__ == '__main__':

    main()

