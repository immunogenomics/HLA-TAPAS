#-*- coding: utf-8 -*-

import argparse, textwrap
from MyManhattan.manhattan import HATK_Manhattan



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='MyManhattan',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        The Manhattan Plot

        Supervisor : WansonChoi
        E-mail : wansonchoi@snu.ac.kr

    #################################################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(ex. *.assoc.logistic).\n\n",
                        nargs='+', required=True)
    # parser.add_argument("--plot-label", "-pl", help="\nPlot Label\n\n", default="")
    parser.add_argument("--out", "-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("--hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)

    parser.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n", default="#FF0000")

    parser.add_argument("--point-size", "-ps"/, help="\nGeneral point size (default: 15).\n\n", default="15")
    parser.add_argument("--yaxis-unit", "-yau", help="\nY-axis value(-log10(x)) unit (default : 10).\n\n", default="10")

    parser.add_argument("--HLA", help="\nWhich HLA gene for Manhattan plot .\n\n",
                        choices=['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')

    ### for Testing
    # args = parser.parse_args([])

    ### for Publish
    args = parser.parse_args()
    print(args)

    HATK_Manhattan(args.assoc_result, args.out, args.hg, _top_color=args.top_color,
                   _point_size=args.point_size, _yaxis_unit=args.yaxis_unit, _HLA=args.HLA)
