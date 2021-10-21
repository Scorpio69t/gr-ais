#!/usr/bin/env python3
# ais decoder, main app

from gnuradio.eng_option import eng_option
from optparse import OptionParser
import ais

import time
import sys
import socket


def main():
    # Create Options Parser
    parser = OptionParser(option_class=eng_option, conflict_handler="resolve")
    ais.radio.ais_radio.add_radio_options(parser)
    (options, args) = parser.parse_args()

    tb = ais.radio.ais_radio(options)
    tb.run()
    tb.close()


if __name__ == '__main__':
    main()
