from IPython.core.magic import register_cell_magic
import string
# import re
import argparse
import pathlib


class NamespaceFormatter(string.Formatter):
   def __init__(self, namespace={}):
       string.Formatter.__init__(self)
       self.namespace = namespace

   def get_value(self, key, args, kwds):
       if isinstance(key, str):
           try:
               # Check explicitly passed arguments first
               return kwds[key]
           except KeyError:
               return self.namespace[key]
       else:
           string.Formatter.get_value(key, args, kwds)


@register_cell_magic
def genfile(line, cell):
    """ Generate file 
    
    This magic is an extension of the %%writefile magic. The additional functionality 
    includes:
      - variable expansion (similar to string.format and f-strings)
      - set file permisssions with -m | --mode flag"""

    # parse the command line arguments
    parser = argparse.ArgumentParser(description='Generate file magic')
    parser.add_argument('outfile', type=str, help='output file')
    parser.add_argument('-m', '--mode', help='file permissions')
    args = parser.parse_args(line.split())

    # format the content
    fmt = NamespaceFormatter(globals())
    content = fmt.format(cell)
    content = content.replace("\n    ", "\n\t")

    # TODO: this should check if the file content changed to
    # prevent updating the timestamp if nothing changed.
    # This will help with Makefile dependencies.

    # generate the file
    filename = pathlib.Path(args.outfile)
    with open(filename, "w") as out:
        out.write(content)
    if args.mode is not None:
        filename.chmod(int(args.mode, 8))

    print(f"Generated file '{filename}'")

    return
