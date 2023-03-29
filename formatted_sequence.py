from matplotlib import colors
import seaborn as sns


class FormattedSequence:

    def __init__(self, seq):
        self.seq = seq
        self._rebuild_array()
        self.tagList = []
        self.formatClassDict = {}
        
        self.palette = sns.color_palette('deep')
        alpha = 0.5
        backgroundColor0 = colors.to_hex(self.palette[0] + (alpha,), keep_alpha=True)
        backgroundColor1 = colors.to_hex(self.palette[1] + (alpha,), keep_alpha=True)
        backgroundColor2 = colors.to_hex(self.palette[2] + (alpha,), keep_alpha=True)
        backgroundColorGrey = colors.to_hex(tuple(3*[0.5]) + (alpha,), keep_alpha=True)
        
        # Predefined format classes
        name = 'start_codon'
        cssAttributes = {'background-color':backgroundColor0, 'font-weight':'bold'}
        self.formatClassDict[name] = cssAttributes
        
        name = 'stop_codon'
        cssAttributes = {'background-color':backgroundColor1}
        self.formatClassDict[name] = cssAttributes
        
        name = 'ShineDalgarno'
        cssAttributes = {'text-decoration':'underline'}
        self.formatClassDict[name] = cssAttributes

        
    def _rebuild_array(self):
        self.arr = [{'pos':i, 'character':c, 'tagList':{'class':[]}} for i, c in enumerate(list(self.seq))]
        
        
    def add_tag_list(self, tagList):
        # Apply tags to the array
        for tag in tagList:
            for pos in range(tag['start'], tag['end']):
                self.arr[pos]['tagList']['class'].append(tag['class'])

        if type(tagList) is list:
            self.tagList += tagList
        elif type(tagList) is dict:
            self.tagList += [tagList]
        
        
    def generate_html_formatted_sequence(self):
        s = ''
        s += '<div class="prettySeq">'
        for d in self.arr:
            tag = 'span'
            tagOpen = '<' + tag + ' class="{}">'.format(' '.join(d['tagList']['class']))
            tagClose = '</' + tag + '>'
            if len(d['tagList']['class']) > 0:
                s += tagOpen + d['character'] + tagClose
            else:
                s += d['character']
        s += '</div>'
        return s
    
    
    def header(self):
        formatClassCssStringDict = {}
        for k, cssAttDict in self.formatClassDict.items():
            cssAtt = ' '.join(['{}: {};'.format(k, v) for k, v in cssAttDict.items()])
            formatClassCssStringDict[k] = cssAtt
        print(formatClassCssStringDict)
        
        style = ("<style>\n    div.prettySeq { font-family:monospace;  word-break: break-all; }\n" +
                 '\n'.join(['    .{} {{ {} }}'.format(name, cssAtt) for name, cssAtt in formatClassCssStringDict.items()]) +
                 "\n</style>\n")
        header = style
        return header
        
        
    def generate_html(self):
        s = ''
        s += '<html>'
        s += self.header()
        s += self.generate_html_formatted_sequence()
        s += '</html>'
        return s