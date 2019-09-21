#   This is a component of LinuxCNC
#   Copyright 2011, 2013 Dewey Garrett <dgarrett@panix.com>, Michael
#   Haberler <git@mah.priv.at>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
import os
import re

class EmcToolTable(object):
    ''' intended as bug-compatible Python replacement for the
    tooltable io used in iocontrol
    NB: old file formats not supported.
    '''

    ttype = { 'T' : int, 'P': int, 'Q':int,
              'X' : float, 'Y' : float, 'Z' : float,
              'A' : float, 'B' : float, 'C' : float,
              'U' : float, 'V' : float, 'W' : float,
              'I' : float, 'J' : float, 'D' : float }

    def __init__(self,filename,random_toolchanger):
         self.filename = filename
         self.random_toolchanger = random_toolchanger

    def load_table(self, tooltable,comments,fms):
        self.fakepocket = 0
        fp = open(self.filename)
        lno = 0
        for line in fp.readlines():
            lno += 1
            if not line.startswith(';') and line.strip():
                entry = self.parseline(lno,line.strip())
                if entry:
                    self.assign(tooltable,entry,comments,fms)
        fp.close()

    def save_table(self, tooltable, comments,fms):
        os.rename(self.filename,self.filename + '.bak')
        fp = open(self.filename, 'w')
        start = 0 if self.random_toolchanger else 1
        for p in range(start,len(tooltable)):
            t = tooltable[p]
            if t.toolno != -1:
                print("T%d P%d" % (t.toolno, p if self.random_toolchanger else fms[p]), end=' ', file=fp)
                if t.diameter:  print("D%f" % (t.diameter), end=' ', file=fp)
                if t.offset.x: print("X%+f" % (t.offset.x), end=' ', file=fp)
                if t.offset.y: print("Y%+f" % (t.offset.y), end=' ', file=fp)
                if t.offset.z: print("Z%+f" % (t.offset.z), end=' ', file=fp)
                if t.offset.a: print("A%+f" % (t.offset.a), end=' ', file=fp)
                if t.offset.b: print("B%+f" % (t.offset.b), end=' ', file=fp)
                if t.offset.c: print("C%+f" % (t.offset.c), end=' ', file=fp)
                if t.offset.u: print("U%+f" % (t.offset.u), end=' ', file=fp)
                if t.offset.v: print("V%+f" % (t.offset.v), end=' ', file=fp)
                if t.offset.w: print("W%+f" % (t.offset.w), end=' ', file=fp)
                if t.offset.w: print("W%+f" % (t.offset.w), end=' ', file=fp)
                if t.frontangle: print("I%+f" % (t.frontangle), end=' ', file=fp)
                if t.backangle: print("J%+f" % (t.backangle), end=' ', file=fp)
                if t.orientation: print("Q%+d" % (t.orientation), end=' ', file=fp)
                if p in comments and comments[p]:
                    print(";%s" % (comments[p]), file=fp)
                else:
                    print(file=fp)
        fp.close()

    def assign(self,tooltable,entry,comments,fms):
        pocket = entry['P']
        if not self.random_toolchanger:
            self.fakepocket += 1
            if self.fakepocket >= len(tooltable):
                print("too many tools. skipping tool %d" % (toolno))
                return
            if not fms is None:
                fms[self.fakepocket] = pocket
            pocket = self.fakepocket
        if pocket < 0 or pocket > len(tooltable):
            print("max pocket number is %d. skipping tool %d" % (len(tooltable) - 1, toolno))
            return

        tooltable[pocket].zero()
        for (key,value) in list(entry.items()):
            if key == 'T' : tooltable[pocket].toolno = value
            if key == 'Q' : tooltable[pocket].orientation = value
            if key == 'D' : tooltable[pocket].diameter = value
            if key == 'I' : tooltable[pocket].frontangle = value
            if key == 'J' : tooltable[pocket].backangle = value
            if key == 'X' : tooltable[pocket].offset.x = value
            if key == 'Y' : tooltable[pocket].offset.y = value
            if key == 'Z' : tooltable[pocket].offset.z = value
            if key == 'A' : tooltable[pocket].offset.a = value
            if key == 'B' : tooltable[pocket].offset.b = value
            if key == 'C' : tooltable[pocket].offset.c = value
            if key == 'U' : tooltable[pocket].offset.u = value
            if key == 'V' : tooltable[pocket].offset.v = value
            if key == 'W' : tooltable[pocket].offset.w = value
            if key == 'comment' : comments[pocket] = value # aaargh

    def parseline(self,lineno,line):
        """
        read a tooltable line
        if an entry was parsed successfully, return a  Tool() instance
        """
        line.rstrip("\n")
        if re.match('\A\s*T\d+',line):
            semi = line.find(";")
            if semi != -1:
                comment = line[semi+1:]
            else:
                comment = None
            entry = line.split(';')[0]
            result = dict()
            for field in entry.split():
                (name,value)  = re.search('([a-zA-Z])([+-]?\d*\.?\d*)',field).groups()
                if name:
                    key = name.upper()
                    result[key] = EmcToolTable.ttype[key](value)
                else:
                    print("%s:%d  bad line: '%s' " % (self.filename, lineno, entry))
            result['comment'] = comment
            return result
        print("%s:%d: unrecognized tool table entry   '%s'" % (self.filename,lineno,line))


    def restore_state(self,e):
        pass

    def save_state(self,e):
        pass
