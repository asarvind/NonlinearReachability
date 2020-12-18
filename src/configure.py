import tkinter as tk
from tkinter import font as tkFont
import tkinter.scrolledtext as st

def configure( mywin, config ):
    refPos = [ 10, 10 ]

    def onoffFn( state ):
        config[ "state" ] = state
        
    onoffLabel = tk.Label( mywin, text = "use server", bg = "white",
                           font = tkFont.Font( size = 22 )
                           )
    onoffLabel.place( x = refPos[0]+80, y = refPos[1]+0 )
    onoff = { "yes", "no" }
    opt = tk.StringVar( mywin )
    opt.set( config[ "state" ] )
    onoffMenu = tk.OptionMenu( mywin, opt, "no", "yes",
                               command = onoffFn )
    onoffMenu[ "font" ] = tkFont.Font( size = 25 )
    onoffMenu[ "bg" ] = "white"
    onoffMenu.place( x = refPos[0], y = refPos[1] )

    sufEntry = tk.Text( mywin, bg = "white", highlightbackground = "black",
                         font = tkFont.Font( size = 17 ), height = 1, width = 30,
                          )
    sufEntry.delete( "1.0", tk.END )
    sufEntry.insert( tk.INSERT, config[ "suffix" ] )
    sufEntry.place( x = refPos[0], y = refPos[1] + 60 )
    sufLabel = tk.Label( mywin, text = "ssh shortcut name", bg = "white",
                           font = tkFont.Font( size = 22 ) )
    sufLabel.place( x = refPos[0], y = refPos[1]+90 )

    pathEntry = tk.Text( mywin, bg = "white", highlightbackground = "black",
                         font = tkFont.Font( size = 17 ), height = 1, width = 60,
                          )
    pathEntry.delete( "1.0", tk.END )
    pathEntry.insert( tk.INSERT, config[ "path" ] )
    pathEntry.place( x = refPos[0], y = refPos[1] + 160 )
    pathLabel = tk.Label( mywin, text = "host path to folder containing tool library", bg = "white",
                           font = tkFont.Font( size = 22 ) )
    pathLabel.place( x = refPos[0], y = refPos[1]+190 )

    

    def save():        
        config[ "path" ] = pathEntry.get( "1.0", tk.END )
        config[ "path" ] = config[ "path" ].replace( " ", "" )
        config[ "path" ] = config[ "path" ].replace( "\n", "" )

        config[ "suffix" ] = sufEntry.get( "1.0", tk.END )
        config[ "suffix" ] = config[ "suffix" ].replace( " ", "" )
        config[ "suffix" ] = config[ "suffix" ].replace( "\n", "" )

        mywin.destroy()

    saveButton =   tk.Button( mywin, text = "save", command = save,
                      font = tkFont.Font( size = 20 ),
                      highlightbackground = "green",
                      bg = "white")        
    saveButton.place( x = refPos[0] + 320, y = refPos[1] + 250 )
    
