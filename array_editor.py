# array_editor.py -- Example of using array editors
import numpy as np
from enthought.traits.api import HasPrivateTraits, Array
from enthought.traits.ui.api import View, ArrayEditor, Item
from enthought.traits.ui.menu import NoButtons

class ArrayEditorTest(HasPrivateTraits):
    three = Array(np.int, (3,3))
    four = Array(np.float,
                 (4,4),
                 editor = ArrayEditor(width = -50))
    view = View(
        Item('three',label='3x3 Integer'),
        '_',
        Item('three',label='Integer Read-only',style='readonly'),
        '_',
        Item('four', label='4x4 Float'),
        '_',
        Item('four', label='Float Read-only',style='readonly'),
        buttons = NoButtons,
        resizable = True )
    
if __name__ == '__main__': ArrayEditorTest().configure_traits()
