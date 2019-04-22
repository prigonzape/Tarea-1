import argparse

def getParser():
    parser = argparse.ArgumentParser(description='Descripcion de la funcion de este programa.')
    parser.add_argument('-o',type=str,dest="formato",help="'png')
    if len(sys.argv) == 1:
        print >> sys.stderr,parser.print_help()
        exit(0)
    return parser

def main():
    args=getParser().parse_args()
    formato=args.formato

  
    if formato=="png":
        plt.savefig("ajustar.png")
    elif formato=="pdf":
        plt.savefig("ajustar.pdf")

