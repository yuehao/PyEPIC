import argparse
import ply.lex as lex
import ply.yacc as yacc
import copy
import re
import numpy



 
class EpicLexer(object):
    literals = "=,+-*/:&()"
    reserved = {'beam1':'KWB1','beam2':'KWB2','end':'KWEND', 'global':'KWGLOBAL'}
    tokens = ['NUMBER' , 'STRING', 'RET']+list(reserved.values())
    
    def t_NUMBER(self, t):
        #r'([+-]?[0-9]*\.?[0-9]+|[0-9]+\.[0-9]*)(([eE][-+]?[0-9]+)?)'
        r'([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
        t.value = float(t.value)
        return t

    t_ignore = ' \t'

    def t_STRING(self, t):
        r'[A-Za-z][A-Za-z0-9_\.]*'
        t.type = EpicLexer.reserved.get(t.value, 'STRING')
        return t

    def t_QUOTE(self, t):
        r'\"\''
        pass

    def t_COMMENT(self, t):
        r'\#.*'
        pass
    
    def t_RET(self, t):
        r'\n+'
        t.lexer.lineno += len(t.value)
        t.value = 'RET'
        return t

    def t_EOF (self, t):
        r'<<EOF>>'
        t.value = 'EOF'
        return t

    def t_error(self, t):
        print("Illegal character {}".format(t.value[0]));
        t.lexer.skip(1)
    
    def __init__(self,**kwargs):
        self.lexer = lex.lex(module=self, reflags=re.IGNORECASE, **kwargs)
        self.continuetor = 0

    def lexit(self ,data):
        self.lexer.input(data)
        while True:
            tok = self.lexer.token()
            if not tok: break
            #print tok


class EpicParser(object):
    #starting = "statement"
    tokens=EpicLexer.tokens

    def __init__(self, **kwargs):
        self.which_beam = 0
        self.temp_dict = {}
        self.nlist = []
        self.parser = yacc.yacc(module=self, debug=0, **kwargs)
    
    def p_statment(self, p):
        '''statment : '&' KWB1 RET
                    | '&' KWB2 RET
                    | '&' KWEND RET     
                    | '&' KWGLOBAL RET
                    | RET
                    | STRING '=' numlist RET
                    | STRING '=' STRING RET
        '''
        if len(p) == 4:
            if p[2] == 'beam1':
                self.which_beam = 1
            elif p[2] == 'beam2':
                self.which_beam = 2
            elif p[2] == 'global':
                self.which_beam = -1
            else:
                self.which_beam = 0
        if len(p) == 5:
            if type(p[3]) is str:
                self.temp_dict[p[1]] = p[3]
            else:
                if len(self.nlist) > 1:
                    self.temp_dict[p[1]] = numpy.array(copy.deepcopy(self.nlist))
                else:
                    self.temp_dict[p[1]] = self.nlist[0]
                self.nlist = []
    
    def p_numlist(self, p):
        '''numlist : numbers
                   | numlist ',' numbers'''
        if len(p) == 2:
            self.nlist = [p[1],]
        if len(p) == 4:
            self.nlist.append(p[3])

    def p_numbers(self, p):
        '''numbers : NUMBER
                   | NUMBER '+' NUMBER
                   | NUMBER '-' NUMBER
                   | NUMBER '*' NUMBER
                   | NUMBER '/' NUMBER'''
        if len(p) == 2:
            p[0] = p[1]
        if len(p) == 4:
            if p[2] == '+':
                p[0] = p[1] + p[3]
            elif p[2] == '-':
                p[0] = p[1] - p[3]
            elif p[2] == '*':
                p[0] = p[1] * p[3]
            elif p[2] == '/':
                p[0] = p[1] / p[3]


    def p_error(self, p):
        raise Exception(p)               


def parse_file(file_stream):
    beam1_dict = {}
    beam2_dict = {}
    global_dict = {}
    which_beam = 0
    for line in file_stream:
        m = EpicLexer()
        m.lexit(line)
        p = EpicParser()
        p.parser.parse(line,lexer=m.lexer)
        if p.which_beam is not 0:
            which_beam = p.which_beam
        if which_beam == 0:
            print("ERROR in the file, no beam is specified")
            exit(-1)
        if which_beam == 1:
            beam1_dict.update(p.temp_dict)
        elif which_beam == 2:
            beam2_dict.update(p.temp_dict)
        else:
            global_dict.update(p.temp_dict)

    return global_dict,beam1_dict,beam2_dict
            
            
            
if __name__ == '__main__':
    parse_file('testinput_example.in')
else:
    epic_parser = argparse.ArgumentParser()
    epic_parser.add_argument('-i', '--input_file', help="input file of beam parameter")
    epic_parser.add_argument('-sw', '--strong_weak', help='treat beam 1 strong', action='store_true', default=True)
    
