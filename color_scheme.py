import matplotlib.colors as mcolors
from Bio import Alphabet


class ColorScheme(object):
    
    def __init__(self, 
                 groups = [], 
                 title = "", 
                 description = "",
                 default_color = [0, 0, 0], 
                 alphabet = None) :
        """  """
        self.title= title
        self.description = description
        self.default_color = default_color
        self.groups = groups
        if alphabet is None:
            alphabet = Alphabet.SingleLetterAlphabet()
            alphabet.letters = ''.join([letter for cg in groups for letter in cg.symbols])
        self.alphabet = alphabet
            
        color = {}

        if alphabet.size > 1:
            # multi-letter symbols, (codons, etc)
            for cg in groups:
                color[cg.symbols] = cg.color
                if cg.symbols not in alphabet.letters:
                    raise KeyError("Colored symbol does not exist in alphabet.")
        else:
            # 1-letter symbols (amino acid, nucleotide)
            for cg in groups :
                for s in cg.symbols :
                    color[s] = cg.color
                    if s not in alphabet.letters:
                        raise KeyError("Colored symbol does not exist in alphabet.")

        self._color = color

    def color(self, symbol) :
        if symbol in self._color :
            return self._color[symbol]
        return self.default_color


        
class ColorGroup(object) :
    """Associate a group of symbols with a color"""
    def __init__(self, symbols, color, description=None) :
        self.symbols = symbols
        if type(color) is list:
            if color[0] > 1 or color[1] > 1 or color[2] > 1:
                color = [c/256 for c in color]
        self.color = mcolors.to_rgb(color)
        self.description = description

         

# From makelogo
nucleotide = ColorScheme([
    ColorGroup("G", "orange"),
    ColorGroup("TU", "red"),
    ColorGroup("C",  "blue"),
    ColorGroup("A",  "green") ])

base_pairing = ColorScheme([
    ColorGroup("TAU",  "darkorange", "Weak (2 Watson-Crick hydrogen bonds)"),
    ColorGroup("GC",    "blue", "Strong (3 Watson-Crick hydrogen bonds)") ])

# From Crooks2004c-Proteins-SeqStr.pdf
hydrophobicity = ColorScheme([
    ColorGroup( "RKDENQ",   "blue", "hydrophilic"),
    ColorGroup( "SGHTAP",   "green", "neutral"  ),
    ColorGroup( "YVMCLFIW", "black",  "hydrophobic") ])

# from makelogo
chemistry = ColorScheme([
  ColorGroup( "GSTYC",  "green",   "polar"),
  ColorGroup( "NQ",      "purple", "neutral"), 
  ColorGroup( "KRH",     "blue",   "basic"),
  ColorGroup( "DE",      "red",    "acidic"),
  ColorGroup("PAWFLIMV", "black",  "hydrophobic") ]
  )   

charge = ColorScheme([
    ColorGroup("KRH", "blue", "Positive" ),
    ColorGroup( "DE", "red", "Negative") ]
    )


taylor = ColorScheme([
    ColorGroup( 'A', '#CCFF00' ),
    ColorGroup( 'C', '#FFFF00' ),
    ColorGroup( 'D', '#FF0000'),
    ColorGroup( 'E', '#FF0066' ),
    ColorGroup( 'F', '#00FF66'),
    ColorGroup( 'G', '#FF9900'),
    ColorGroup( 'H', '#0066FF'),
    ColorGroup( 'I', '#66FF00'),
    ColorGroup( 'K', '#6600FF'),
    ColorGroup( 'L', '#33FF00'),
    ColorGroup( 'M', '#00FF00'),
    ColorGroup( 'N', '#CC00FF'),
    ColorGroup( 'P', '#FFCC00'),
    ColorGroup( 'Q', '#FF00CC'),
    ColorGroup( 'R', '#0000FF'),
    ColorGroup( 'S', '#FF3300'),
    ColorGroup( 'T', '#FF6600'),
    ColorGroup( 'V', '#99FF00'),
    ColorGroup( 'W', '#00CCFF'),
    ColorGroup( 'Y', '#00FFCC')],
    title = "Taylor",
    description = "W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997)"
    )

rasmol = ColorScheme([
    ColorGroup('A', [200, 200, 200]),
    ColorGroup('R', [20, 90, 255]),
    ColorGroup('N', [0, 220, 220]),
    ColorGroup('D', [230, 10, 10]),
    ColorGroup('C', [230, 230, 0]),
    ColorGroup('Q', [0, 220, 220]),
    ColorGroup('E', [230, 10, 10]),
    ColorGroup('G', [235, 235, 235]),
    ColorGroup('H', [130, 130, 210]),
    ColorGroup('I', [15, 130, 15]),
    ColorGroup('L', [15, 130, 15]),
    ColorGroup('K', [20, 90, 255]),
    ColorGroup('M', [230, 230, 0]),
    ColorGroup('F', [50, 50, 170]),
    ColorGroup('P', [220, 150, 130]),
    ColorGroup('S', [250, 150, 0]),
    ColorGroup('T', [250, 150, 0]),
    ColorGroup('W', [180, 90, 180]),
    ColorGroup('Y', [50, 50, 170]),
    ColorGroup('V', [15, 130, 15]),
    ColorGroup('B', [255, 105, 180]),
    ColorGroup('Z', [255, 105, 180]),
    ColorGroup('X', [190, 160, 110])
    ],
    default_color='none',
    title = "rasmol",
    description = ""
    )