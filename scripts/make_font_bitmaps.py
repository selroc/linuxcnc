#!/usr/bin/python2

# 1) Make bitmap files for font for use with python3 which
#    lacks some pango/cairo/pangocairo capabilities
# 2) Adapted from glnav.py:use_pango_font()
# 3) Requires python2 for availability of pangocairo
# 4) Example failure with python3 using gi.repository:
#    AttributeError: 'gi.repository.PangoCairo' object
#                     has no attribute 'CairoContext'

import sys
import array
import itertools
import pango
import cairo
import pangocairo

def write_bitmap(fontname,start,count,DIM,tofile):
    dmax=wmax=hmax=0
    fontDesc = pango.FontDescription(fontname)
    # fdata is the font bitmap data b:signed, B:unsigned
    # NOTE: 'B'
    fdata = array.array('B', itertools.repeat(0, DIM*DIM))
    surface = cairo.ImageSurface.create_for_data(fdata, cairo.FORMAT_A8, DIM, DIM)
    context = pangocairo.CairoContext(cairo.Context(surface))
    layout = context.create_layout()
    fontmap = pangocairo.cairo_font_map_get_default()
    font = fontmap.load_font(fontmap.create_context(), fontDesc)
    layout.set_font_description(fontDesc)
    metrics = font.get_metrics()
    descent = metrics.get_descent()
    d = pango.PIXELS(descent)
    linespace = metrics.get_ascent() + metrics.get_descent()
    width = metrics.get_approximate_char_width()

    f = open(tofile,'wb')
    f.write("""# Created by: %s
# fontname  = '%s'
# start     = %d
# count     = %d
# DIM       = %d
# tofile    = %s

# NOTE: Edit this file to remove unneeded fbits dictionary
#       entries for unneeded characters

fbits=dict()
"""%(sys.argv[0],fontname,start,count,DIM,tofile))
    for i in range(count):
        ch = unichr(start+i)
        layout.set_text(ch)
        w, h = layout.get_size()
        context.save()
        context.new_path()
        context.rectangle(0, 0, DIM, DIM)
        context.set_source_rgba(0., 0., 0., 0.)
        context.set_operator (cairo.OPERATOR_SOURCE);
        context.paint()
        context.restore()

        context.save()
        context.set_source_rgba(1., 1., 1., 1.)
        context.set_operator (cairo.OPERATOR_SOURCE);
        context.move_to(0, 0)
        context.update_layout(layout)
        context.show_layout(layout)
        context.restore()

        if i>=32 and i<=90: # punctuation,numbers,uppercase
            w, h = pango.PIXELS(w), pango.PIXELS(h)
            if (w > wmax): wmax=w;wmaxchar=ch
            if (h > hmax): hmax=h;hmaxchar=ch

        f.write("#index=%2d, char=<%s>\n"%(i,str(ch)))
        f.write("fbits[%d]=B'\\\n"%i)
        for row in range(0,DIM):
            for col in range(0,DIM):
                v=fdata[row*DIM +col]
                v='%02x'%v
                f.write("\\x%s"%v)
            f.write("\\\n")
        f.write("'\n")
    f.write("""
width   = %2d
height  = %2d
descent = %2d
DIM     = %2d
"""%(wmax,hmax,dmax,DIM))
    #print("wmax=%2d"%wmax)
    #print("hmax=%2d"%hmax)
    f.close()
#---------------------------------------------------------------

# edit as required:
fontname  = "courier bold 11"
start     =   0
count     = 128
DIM       =  32
tofile    = '/tmp/%s.bitmap'%fontname.replace(' ','_')

print("fontname =%s"%fontname)
print("start    =%3d"%start)
print("count    =%3d"%count)
print("DIM      =%3d"%DIM)
print("tofile   =%s"%tofile)

write_bitmap(fontname,start,count,DIM,tofile)
