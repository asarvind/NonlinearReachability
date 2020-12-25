import pickle
import os

from configure import *
from symprocess import *


#----------------------------------------------------------------------
# create UI
#----------------------------------------------------------------------
os.system( "mkdir -p src/pywrite" )
os.system( "mkdir -p results" )
ui = tk.Tk( )
ui.title( "Flowpipe Simulator" )
ui.geometry( "500x220" )
ui.configure( bg="white" )

header = tk.Label( ui, text = "Flowpipe Simulator",
                   fg = "green", bg = "white" )
#header.place( x = 500, y = 5 )
header[ "font" ] = tkFont.Font( size = 30, family = "zapfino" )

startPos = [10,20]
yGap = 50
#----------------------------------------------------------------------
# declare list of variables and dynamics (dictionary)
#----------------------------------------------------------------------
global stBounds 
global inpBounds 
global parEqs 
global dynamics

stBounds = {}
inpBounds = {}
parEqs = {}
dynamics = {}



#----------------------------------------------------------------------
# edit state
#----------------------------------------------------------------------
class stWin:
    def __init__(self):
        refPos = [10,20]
        self.ui = tk.Toplevel( ui )
        self.ui.title( "Edit State" )
        self.ui.geometry( "400x650" )
        self.ui.configure( bg = "white" )

        #----------------------------------------------------------------------
        # buttons
        #----------------------------------------------------------------------
        self.closeButton = tk.Button( self.ui, text = "ok", command = self.close,
                                 highlightbackground = "green",
                                 bg = "white")
        self.closeButton[ "font" ] = tkFont.Font( size = 15 )
        self.closeButton.place( x = 200, y = 610 )
        self.stCaption = tk.Label( self.ui, text = "state\n variable\n", bg = "white" )
        self.stCaption[ "font" ] = tkFont.Font( size = 20 )
        self.stCaption.place( x=refPos[0], y=refPos[1]-15 )

        self.stEntry = tk.Text( self.ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 6,
                           highlightbackground = "black", bg = "white")
        self.stEntry.place( x=refPos[0]+10, y=refPos[1]+40 )

        self.stLbCaption = tk.Label( self.ui, text = "initial lower\n bound\n", bg = "white" )
        self.stLbCaption[ "font" ] = tkFont.Font( size = 20 )
        self.stLbCaption.place( x=refPos[0]+110, y=refPos[1]-15 )

        self.stLbEntry = tk.Text( self.ui, font = tkFont.Font( size = 15 ),
                           height = 1, width = 10,
                           highlightbackground = "black", bg = "white")
        self.stLbEntry.place( x=refPos[0]+110, y=refPos[1]+40 )

        self.stUbCaption = tk.Label( self.ui, text = "initial upper\n bound\n", bg = "white" )
        self.stUbCaption[ "font" ] = tkFont.Font( size = 20 )
        self.stUbCaption.place( x=refPos[0]+230, y=refPos[1]-15 )

        self.stUbEntry = tk.Text( self.ui, font = tkFont.Font( size = 15 ),
                           height = 1, width = 10,
                           highlightbackground = "black", bg = "white")
        self.stUbEntry.place( x=refPos[0]+240, y=refPos[1]+40 )

        self.stDisplay = st.ScrolledText( self.ui, font = tkFont.Font( size = 17 ),
                                     height = 22, width = 32,
                                     highlightbackground = "black", bg = "white")
        self.stDisplay.place( x=refPos[0]+15, y=refPos[1]+118 )
        self.stDisplay.configure( state = "disabled" )
        
        self.stEntryButton = tk.Button( self.ui, text = "add/edit",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.stEntryButtonFn )
        self.stEntryButton.place( x = refPos[0]+35, y = refPos[1]+80  )

        self.stDelButton = tk.Button( self.ui, text = "del",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.stDelButtonFn )
        self.stDelButton.place( x = refPos[0]+140, y = refPos[1]+80  )

        self.stClButton = tk.Button( self.ui, text = "clear",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.stClButtonFn )
        self.stClButton.place( x = refPos[0]+220, y = refPos[1]+80 )

        # load state
        self.stDisplay.configure( state = "normal" )
        for x in stBounds:
            self.stDisplay.insert( tk.INSERT, x + " : " + str(stBounds[ x ])+"\n" )
        self.stDisplay.configure( state = "disabled" )


    def close(self):
        self.ui.destroy()

    # insert button
    def stEntryButtonFn(self):
        v = self.stEntry.get( "1.0", tk.END )
        v = v.replace( "\n", "" )
        v = v.replace( " ", "" )

        lb = self.stLbEntry.get( "1.0", tk.END )
        lb = lb.replace( "\n", "" )
        lb = lb.replace( " ", "" )

        ub = self.stUbEntry.get( "1.0", tk.END )
        ub = ub.replace( "\n", "" )
        ub = ub.replace( " ", "" )

        if lb == "" and ub != "":
            lb = ub
        elif ub == "" and lb != "":
            ub = lb
        elif ub == "" and lb == "":
            lb = str( 0 )
            ub = str( 0 )

        if v != "":
            stBounds[ v ] = [ float( lb ), float( ub ) ]
            self.stDisplay.configure( state = "normal" )
            self.stDisplay.delete( "1.0", tk.END )
            # stEntry.delete( "1.0", tk.END )
            # stLbEntry.delete( "1.0", tk.END )
            # stUbEntry.delete( "1.0", tk.END )
            for x in stBounds:
                self.stDisplay.insert( tk.INSERT, x + " : " + str(stBounds[ x ])+"\n" )
            self.stDisplay.configure( state = "disabled" )
        
    # delete button
    def stDelButtonFn(self):
        v = self.stEntry.get( "1.0", tk.END )
        v = v.replace( "\n", "" )
        v = v.replace( " ", "" )
        if v != "" and ( v in stBounds.keys() ):
            stBounds.pop( v )
            self.stDisplay.configure( state = "normal" )
            self.stDisplay.delete( "1.0", tk.END )
            self.stEntry.delete( "1.0", tk.END )
            self.stLbEntry.delete( "1.0", tk.END )
            self.stUbEntry.delete( "1.0", tk.END )
            for x in stBounds:
                self.stDisplay.insert( tk.INSERT, x + " : " + str(stBounds[ x ])+"\n" )
            self.stDisplay.configure( state = "disabled" )

    # clear button
    def stClButtonFn(self):
        self.stDisplay.configure( state = "normal" )
        stBounds.clear()
        self.stDisplay.delete( "1.0", tk.END )
        self.stDisplay.configure( state = "disabled" )
        self.stEntry.delete( "1.0", tk.END )
        self.stLbEntry.delete( "1.0", tk.END )
        self.stUbEntry.delete( "1.0", tk.END )
            



editStButton = tk.Button( ui, text = "State", command = stWin,
                          font = tkFont.Font( size = 25 ),
                          highlightbackground = "yellow",
                          bg = "white")
editStButton.place( x = startPos[0], y = startPos[1]  )


class editstate( stWin ):
    def __init__( self ):
        stWin.__init__( self )
        self.stDelButton.destroy()
        self.stClButton.destroy()
    
#----------------------------------------------------------------------
# edit input
#----------------------------------------------------------------------
class inpWin:
    def __init__(self):
        refPos = [10,20]
        self.ui = tk.Toplevel( ui )
        self.ui.title( "Edit Input" )
        self.ui.geometry( "400x650" )
        self.ui.configure( bg = "white" )

        #----------------------------------------------------------------------
        # buttons
        #----------------------------------------------------------------------
        self.closeButton = tk.Button( self.ui, text = "ok", command = self.close,
                                 highlightbackground = "green",
                                 bg = "white")
        self.closeButton[ "font" ] = tkFont.Font( size = 15 )
        self.closeButton.place( x = 200, y = 610 )
        self.inpCaption = tk.Label( self.ui, text = "input\n variable\n", bg = "white" )
        self.inpCaption[ "font" ] = tkFont.Font( size = 20 )
        self.inpCaption.place( x=refPos[0], y=refPos[1]-15 )

        self.inpEntry = tk.Text( self.ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 6,
                           highlightbackground = "black", bg = "white")
        self.inpEntry.place( x=refPos[0]+10, y=refPos[1]+40 )

        self.inpLbCaption = tk.Label( self.ui, text = "lower\n bound\n", bg = "white" )
        self.inpLbCaption[ "font" ] = tkFont.Font( size = 20 )
        self.inpLbCaption.place( x=refPos[0]+110, y=refPos[1]-15 )

        self.inpLbEntry = tk.Text( self.ui, font = tkFont.Font( size = 15 ),
                           height = 1, width = 10,
                           highlightbackground = "black", bg = "white")
        self.inpLbEntry.place( x=refPos[0]+110, y=refPos[1]+40 )

        self.inpUbCaption = tk.Label( self.ui, text = "upper\n bound\n", bg = "white" )
        self.inpUbCaption[ "font" ] = tkFont.Font( size = 20 )
        self.inpUbCaption.place( x=refPos[0]+230, y=refPos[1]-15 )

        self.inpUbEntry = tk.Text( self.ui, font = tkFont.Font( size = 15 ),
                           height = 1, width = 10,
                           highlightbackground = "black", bg = "white")
        self.inpUbEntry.place( x=refPos[0]+240, y=refPos[1]+40 )

        self.inpDisplay = st.ScrolledText( self.ui, font = tkFont.Font( size = 17 ),
                                     height = 22, width = 32,
                                     highlightbackground = "black", bg = "white")
        self.inpDisplay.place( x=refPos[0]+18, y=refPos[1]+118 )
        self.inpDisplay.configure( state = "disabled" )
        
        self.inpEntryButton = tk.Button( self.ui, text = "add/edit",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.inpEntryButtonFn )
        self.inpEntryButton.place( x = refPos[0]+35, y = refPos[1]+80  )

        self.inpDelButton = tk.Button( self.ui, text = "del",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.inpDelButtonFn )
        self.inpDelButton.place( x = refPos[0]+140, y = refPos[1]+80  )

        self.inpClButton = tk.Button( self.ui, text = "clear",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.inpClButtonFn )
        self.inpClButton.place( x = refPos[0]+220, y = refPos[1]+80 )
        
        # load input
        self.inpDisplay.configure( state = "normal" )
        for u in inpBounds:
            self.inpDisplay.insert( tk.INSERT, u + " : " + str(inpBounds[ u ])+"\n" )
        self.inpDisplay.configure( state = "disabled" )



    def close(self):
        self.ui.destroy()

    # insert button
    def inpEntryButtonFn(self):
        v = self.inpEntry.get( "1.0", tk.END )
        v = v.replace( "\n", "" )
        v = v.replace( " ", "" )

        lb = self.inpLbEntry.get( "1.0", tk.END )
        lb = lb.replace( "\n", "" )
        lb = lb.replace( " ", "" )

        ub = self.inpUbEntry.get( "1.0", tk.END )
        ub = ub.replace( "\n", "" )
        ub = ub.replace( " ", "" )

        if lb == "" and ub != "":
            lb = ub
        elif ub == "" and lb != "":
            ub = lb
        elif ub == "" and lb == "":
            lb = str( 0 )
            ub = str( 0 )

        if v != "":
            inpBounds[ v ] = [ float( lb ), float( ub ) ]
            self.inpDisplay.configure( state = "normal" )
            self.inpDisplay.delete( "1.0", tk.END )
            # self.inpEntry.delete( "1.0", tk.END )
            # self.inpLbEntry.delete( "1.0", tk.END )
            # self.inpUbEntry.delete( "1.0", tk.END )
            for x in inpBounds:
                self.inpDisplay.insert( tk.INSERT, x + " : " + str(inpBounds[ x ])+"\n" )
            self.inpDisplay.configure( state = "disabled" )
        
    # delete button
    def inpDelButtonFn(self):
        v = self.inpEntry.get( "1.0", tk.END )
        v = v.replace( "\n", "" )
        v = v.replace( " ", "" )
        if v != "" and ( v in inpBounds.keys() ):
            inpBounds.pop( v )
            self.inpDisplay.configure( state = "normal" )
            self.inpDisplay.delete( "1.0", tk.END )
            self.inpEntry.delete( "1.0", tk.END )
            self.inpLbEntry.delete( "1.0", tk.END )
            self.inpUbEntry.delete( "1.0", tk.END )
            for x in inpBounds:
                self.inpDisplay.insert( tk.INSERT, x + " : " + str(inpBounds[ x ])+"\n" )
            self.inpDisplay.configure( state = "disabled" )

    # clear button
    def inpClButtonFn(self):
        self.inpDisplay.configure( state = "normal" )
        inpBounds.clear()
        self.inpDisplay.delete( "1.0", tk.END )
        self.inpDisplay.configure( state = "disabled" )
        self.inpEntry.delete( "1.0", tk.END )
        self.inpLbEntry.delete( "1.0", tk.END )
        self.inpUbEntry.delete( "1.0", tk.END )
            

editInpButton = tk.Button( ui, text = "Input", command = inpWin,
                          font = tkFont.Font( size = 23 ),
                          highlightbackground = "yellow",
                          bg = "white")
editInpButton.place( x = startPos[0], y = startPos[1]+yGap  )

class editinp( inpWin ):
    def __init__( self ):
        inpWin.__init__( self )
        self.inpDelButton.destroy()
        self.inpClButton.destroy()

#----------------------------------------------------------------------
# entry field for dynamics
#----------------------------------------------------------------------
class fieldUi:
    def __init__(self):
        refPos = [50,20]
        xPos = refPos[0]
        yPos = refPos[1]

        self.ui = tk.Toplevel( ui )
        self.ui.title( "Edit Dynamics" )
        self.ui.geometry( "800x650" )
        self.ui.configure( bg = "white" )
        
        self.dynamicsDisplay = st.ScrolledText( self.ui, height=16, width = 57,
                                           highlightbackground = "black",
                                           font = tkFont.Font( size = 18 ) )
        self.dynamicsDisplay.place( x = refPos[0], y = yPos+180 )
        self.dynamicsDisplay.configure( state = "disabled" )
        
        # insert state captions
        self.derStVar = tk.Label( self.ui, text="state var\n name", bg="white" )
        self.derStVar.place( x = xPos, y = yPos )
        self.derStVar[ "font" ] = tkFont.Font( size=20 )
        xPos += 100

        self.derCaption = tk.Label( self.ui, text="time derivative expression", bg="white" )
        self.derCaption.place( x=xPos, y=yPos )
        self.derCaption[ "font" ] = tkFont.Font( size=20 )
        xPos += 100

        # reset positions
        xPos = refPos[0]
        yPos += 60

        self.ddt = tk.Label( self.ui, text="d/dt", bg="white",
                        highlightbackground = "white" )
        self.ddt[ "font" ] = tkFont.Font( size=25 )
        self.ddt.place( x=refPos[0]-50, y=refPos[1]+60 )

        self.eqto = tk.Label( self.ui, text="=", bg="white",
                        highlightbackground = "white" )
        self.eqto[ "font" ] = tkFont.Font( size=30 )
        self.eqto.place( x=refPos[0]+80, y=refPos[1]+50 )

        self.derStEntry = tk.Text( self.ui, font=tkFont.Font( size=20 ), height = 1,
                              highlightbackground = "black",  width = 5 )
        self.derStEntry.place( x=xPos, y=yPos )
        xPos += 110

        self.dynamicsEntry = tk.Text( self.ui, font=tkFont.Font( size=18 ), height = 2,
                                  highlightbackground = "black", width = 50 )                          
        self.dynamicsEntry.place( x=xPos, y=yPos-10 )

        self.dynEntryButton = tk.Button( self.ui, text = "add/edit",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 3, font = tkFont.Font( size = 15 ),
                                   command = self.stTable )
        self.dynEntryButton.place( x = refPos[0]+35, y = refPos[1]+125 )

        self.dynDelButton = tk.Button( self.ui, text = "del",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 3, font = tkFont.Font( size = 15 ),
                                   command = self.delDynamics )
        self.dynDelButton.place( x = refPos[0]+140, y = refPos[1]+125  )

        self.dynClButton = tk.Button( self.ui, text = "clear",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 3, font = tkFont.Font( size = 15 ),
                                   command = self.dynClButtonFn )
        self.dynClButton.place( x = refPos[0]+220, y = refPos[1]+125 )

        # initial display
        self.dynamicsDisplay.configure( state = "normal" )
        for ky in dynamics:
            self.dynamicsDisplay.insert( tk.INSERT, "der(" + ky + ") = " + dynamics[ ky ] + "\n" )
        self.dynamicsDisplay.configure( state = "disabled" )

        # ok button
        self.closeButton = tk.Button( self.ui, text = "ok", command = self.close,
                                      highlightbackground = "green",
                                      bg = "white", font = tkFont.Font( size = 20 ))
        self.closeButton.place( x = 400, y = 600 )

    # button to enter dynamics
    def stTable( self ):
        v = self.derStEntry.get( "1.0", "end" )
        v = v.replace("\n","")
        v = v.replace(" ","")
        f = self.dynamicsEntry.get( "1.0", "end" )
        f = f.replace("\n","")
        f = f.replace(" ","")
        if f != "" and v != "":
            dynamics[ v ] = f
            self.dynamicsEntry.delete( "1.0", "end" )
            self.derStEntry.delete( "1.0", "end" )
            self.dynamicsDisplay.configure( state = "normal" )
            self.dynamicsDisplay.delete( "1.0", "end" )
            for ky in dynamics:
                self.dynamicsDisplay.insert( tk.INSERT, "der(" + ky + ") = " + dynamics[ ky ] + "\n" )
            self.dynamicsDisplay.configure( state = "disabled" )

    # button to delete dynamics
    def delDynamics( self ):
        v = self.derStEntry.get( "1.0", "end" )
        v = v.replace("\n","")
        v = v.replace(" ","")
        if v != "":
            dynamics.pop( v )
            self.dynamicsEntry.delete( "1.0", "end" )
            self.derStEntry.delete( "1.0", "end" )
            self.dynamicsDisplay.configure( state = "normal" )
            self.dynamicsDisplay.delete( "1.0", "end" )
            for ky in dynamics:
                self.dynamicsDisplay.insert( tk.INSERT, "der(" + ky + ") = " + dynamics[ ky ] + "\n" )        
            self.dynamicsDisplay.configure( state = "disabled" )

    # clear button
    def dynClButtonFn( self ):
        self.dynamicsDisplay.configure( state = "normal" )
        dynamics.clear()
        self.dynamicsDisplay.delete( "1.0", tk.END )
        self.dynamicsDisplay.configure( state = "disabled" )
        self.dynamicsEntry.delete( "1.0", tk.END )
        self.derStEntry.delete( "1.0", tk.END )

    def close( self ):
        self.ui.destroy()

editDynamicsButton = tk.Button( ui, text = "Dynamics", command = fieldUi,
                          font = tkFont.Font( size = 23 ),
                          highlightbackground = "yellow",
                          bg = "white")
editDynamicsButton.place( x = startPos[0], y = startPos[1]+3*yGap  )

class editdyn( fieldUi ):
    def __init__( self ):
        fieldUi.__init__( self )
        self.dynEntryButton[ "text" ] =  "symbolic processing completed"
        self.dynClButton.destroy()
        self.dynDelButton.destroy()
        
    def stTable( self ):
        pass
        

#----------------------------------------------------------------------
# edit parameter
#----------------------------------------------------------------------
class parWin:
    def __init__(self):
        refPos = [10,20]
        self.ui = tk.Toplevel( ui )
        self.ui.title( "Edit Parameter" )
        self.ui.geometry( "800x650" )
        self.ui.configure( bg = "white" )

        #----------------------------------------------------------------------
        # buttons
        #----------------------------------------------------------------------
        self.closeButton = tk.Button( self.ui, text = "ok", command = self.close,
                                 highlightbackground = "green",
                                 bg = "white")
        self.closeButton[ "font" ] = tkFont.Font( size = 20 )
        self.closeButton.place( x = 400, y = 610 )
        
        self.parCaption = tk.Label( self.ui, text = "parameter\n variable\n", bg = "white" )
        self.parCaption[ "font" ] = tkFont.Font( size = 20 )
        self.parCaption.place( x=refPos[0], y=refPos[1]-15 )

        self.parEntry = tk.Text( self.ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 6,
                           highlightbackground = "black", bg = "white")
        self.parEntry.place( x=refPos[0]+10, y=refPos[1]+40 )

        self.parAsCaption = tk.Label( self.ui, text = "assignment\n", bg = "white" )
        self.parAsCaption[ "font" ] = tkFont.Font( size = 20 )
        self.parAsCaption.place( x=refPos[0]+110, y=refPos[1]-15 )

        self.parAsEntry = tk.Text( self.ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 50,
                           highlightbackground = "black", bg = "white")
        self.parAsEntry.place( x=refPos[0]+110, y=refPos[1]+40 )


        self.parDisplay = st.ScrolledText( self.ui, font = tkFont.Font( size = 20 ),
                                     height = 18, width = 55,
                                     highlightbackground = "black", bg = "white")
        self.parDisplay.place( x=refPos[0]+20, y=refPos[1]+130 )
        self.parDisplay.configure( state = "disabled" )
        
        self.parEntryButton = tk.Button( self.ui, text = "add/edit",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.parEntryButtonFn )
        self.parEntryButton.place( x = refPos[0]+35, y = refPos[1]+80  )

        self.parDelButton = tk.Button( self.ui, text = "del",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.parDelButtonFn )
        self.parDelButton.place( x = refPos[0]+140, y = refPos[1]+80  )

        self.parClButton = tk.Button( self.ui, text = "clear",  fg = "black",
                                   highlightbackground = "green",
                                   bd = 2, font = tkFont.Font( size = 15 ),
                                   command = self.parClButtonFn )
        self.parClButton.place( x = refPos[0]+220, y = refPos[1]+80 )
        
        # load parameter
        self.parDisplay.configure( state = "normal" )
        for u in parEqs:
            self.parDisplay.insert( tk.INSERT, u + " = " + parEqs[ u ] + "\n" )
        self.parDisplay.configure( state = "disabled" )



    def close(self):
        self.ui.destroy()

    # insert button
    def parEntryButtonFn(self):
        v = self.parEntry.get( "1.0", tk.END )
        v = v.replace( "\n", "" )
        v = v.replace( " ", "" )

        eq = self.parAsEntry.get( "1.0", tk.END )
        eq = eq.replace( "\n", "" )
        eq = eq.replace( " ", "" )

        if v != "":
            parEqs[ v ] = eq
            self.parDisplay.configure( state = "normal" )
            self.parDisplay.delete( "1.0", tk.END )
            self.parEntry.delete( "1.0", tk.END )
            self.parAsEntry.delete( "1.0", tk.END )
            for x in parEqs:
                self.parDisplay.insert( tk.INSERT, x + " = " + parEqs[ x ] + "\n" )
            self.parDisplay.configure( state = "disabled" )
        
    # delete button
    def parDelButtonFn(self):
        v = self.parEntry.get( "1.0", tk.END )
        v = v.replace( "\n", "" )
        v = v.replace( " ", "" )
        if v != "" and ( v in parEqs.keys() ):
            parEqs.pop( v )
            self.parDisplay.configure( state = "normal" )
            self.parDisplay.delete( "1.0", tk.END )
            self.parEntry.delete( "1.0", tk.END )
            self.parAsEntry.delete( "1.0", tk.END )
            for x in parEqs:
                self.parDisplay.insert( tk.INSERT, x + " = " + parEqs[ x ] + "\n" )
            self.parDisplay.configure( state = "disabled" )

    # clear button
    def parClButtonFn(self):
        self.parDisplay.configure( state = "normal" )
        parEqs.clear()
        self.parDisplay.delete( "1.0", tk.END )
        self.parDisplay.configure( state = "disabled" )
        self.parEntry.delete( "1.0", tk.END )
        self.parAsEntry.delete( "1.0", tk.END )
            

editParButton = tk.Button( ui, text = "Parameter", command = parWin,
                          font = tkFont.Font( size = 23 ),
                          highlightbackground = "yellow",
                          bg = "white")
editParButton.place( x = startPos[0], y = startPos[1]+2*yGap  )

class editpar( parWin ):
    def __init__( self ):
        parWin.__init__( self )
        self.parClButton.destroy()
        self.parDelButton.destroy()
        self.parEntryButton[ "text" ] =  "symbolic processing completed"

    def parEntryButtonFn( self ):
        pass0


#----------------------------------------------------------------------
# server button
#----------------------------------------------------------------------
server = {}
server[ "state" ] = "no"
server[ "suffix" ] = ""
server[ "path" ] = ""

def serverFn():
    mywin = tk.Toplevel( ui )
    mywin.geometry( "700x350" )
    mywin.configure( bg = "white" )
    configure( mywin, server )
    
serverButton = tk.Button( ui, text = "SSH", command = serverFn,
                          font = tkFont.Font( size = 18 ),
                          highlightbackground = "green",
                          bg = "white" )
serverButton.place( x = startPos[0] + 350, y = startPos[1]  )


#----------------------------------------------------------------------
# Simulate ui
#----------------------------------------------------------------------
maxTimeEntry = tk.Text( ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 5,
                           highlightbackground = "black", bg = "white")
maxTimeEntryLabel = tk.Label( ui, text="max. time", bg="white",
                              font = tkFont.Font( size = 20 ) )


logDivsEntry = tk.Text( ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 5,
                           highlightbackground = "black", bg = "white")
logDivsEntryLabel = tk.Label( ui, text="log_e( no. divisions )", bg="white",
                              font = tkFont.Font( size = 20 ) )

timeStepEntry = tk.Text( ui, font = tkFont.Font( size = 20 ),
                           height = 1, width = 5,
                           highlightbackground = "black", bg = "white" )
timeStepEntryLabel = tk.Label( ui, text="step size", bg="white",
                               font = tkFont.Font( size = 20 ) )

def simulate():
    # write state bounds
    lbobj = open( "src/pywrite/stlb.txt", "w" )
    ubobj = open( "src/pywrite/stub.txt", "w" )
    for mykey in stBounds:
        lbobj.write( str( stBounds[mykey][0] ) + " ")
        ubobj.write( str( stBounds[mykey][1] ) + " ")
    lbobj.close()
    ubobj.close()

    # write input bounds
    lbobj = open( "src/pywrite/inplb.txt", "w" )
    ubobj = open( "src/pywrite/inpub.txt", "w" )
    for mykey in inpBounds:
        lbobj.write( str( inpBounds[mykey][0] ) + " ")
        ubobj.write( str( inpBounds[mykey][1] ) + " ")
    lbobj.close()
    ubobj.close()

    #----------------------------------------------------------------------
    # write simulation parameters
    #----------------------------------------------------------------------
    simwrite = open( "src/pywrite/simpars.txt", "w" )
    
    xstr = maxTimeEntry.get( "1.0", tk.END )
    xstr = xstr.replace( " ", "" )
    xstr = xstr.replace( "\n", "" )
    if xstr != "":
        simstr = xstr + " "
    else:
        simstr = "0.0 "

    xstr = logDivsEntry.get( "1.0", tk.END )
    xstr = xstr.replace( " ", "" )
    xstr = xstr.replace( "\n", "" )
    if xstr != "":
        simstr += xstr + " "
    else:
        simstr = "0.0 "

    xstr = timeStepEntry.get( "1.0", tk.END )
    xstr = xstr.replace( " ", "" )
    xstr = xstr.replace( "\n", "" )
    if xstr != "":
        simstr += xstr + " "
    else:
        simstr = "0.0 "

    simwrite.write( simstr )
    simwrite.close()
    #----------------------------------------------------------------------

    # run executable
    if server[ "state" ] == "no":
        os.system( "make execute" )
    else:
        pth = server[ "path" ]
        execstr = "scp " + "src/pywrite/*.txt " + server[ "suffix" ]
        execstr += ":" + pth + "NonlinearReachability/src/pywrite/"
        os.system( execstr )
        
        # execute inside server
        execstr =  "ssh " + server[ "suffix" ] + " \"cd "
        execstr += pth + "NonlinearReachability/; " + "make execute\"; cd;"
        os.system( execstr )


#----------------------------------------------------------------------
# save and load buttons
#----------------------------------------------------------------------
# entry field for filename
fileEntry = tk.Text( ui, font = tkFont.Font( size = 20 ), bg = "white",
                     highlightbackground = "black", height = 1,
                     width = 15 )
fileEntry.place( x = startPos[0] + 200, y = 2*yGap+30 )
fileExtension = tk.Label( ui, font = tkFont.Font( size = 20 ), bg = "white",
                          text = ".pickle" )
fileExtension.place( x = startPos[0] + 405, y = 2*yGap+30 )
fileEntryLabel = tk.Label( ui, font = tkFont.Font( size = 20 ), bg = "white",
                          text = "filename" )
fileEntryLabel.place( x = startPos[0] + 200, y = 2*yGap-5 )

# save function
def saveInf():
    fname = fileEntry.get( "1.0", tk.END )
    fname = fname.replace( "\n", "" )
    fname = fname.replace( " ", "" )
    fname += ".pickle"
    pick = open( fname, "wb" )
    pickle.dump( [ stBounds, inpBounds, parEqs, dynamics ], pick,
                 protocol=pickle.HIGHEST_PROTOCOL )

# load function
def loadInf():
    fname = fileEntry.get( "1.0", tk.END )
    fname = fname.replace( "\n", "" )
    fname = fname.replace( " ", "" )
    fname += ".pickle"
    pick = open( fname, "rb" )
    global stBounds
    global inpBounds
    global parEqs
    global dynamics
    stBounds, inpBounds, parEqs, dynamics = pickle.load( pick )

saveButton = tk.Button( ui, text = "Save", command = saveInf,
                          font = tkFont.Font( size = 17 ),
                          highlightbackground = "green",
                          bg = "white")
saveButton.place( x = startPos[0] + 200, y = 2*yGap+70 )

loadButton = tk.Button( ui, text = "Load", command = loadInf,
                          font = tkFont.Font( size = 17 ),
                          highlightbackground = "green",
                          bg = "white")
loadButton.place( x = startPos[0] + 300, y = 2*yGap+70 )

#----------------------------------------------------------------------
# process button
#----------------------------------------------------------------------
def process():
    # change geometry of main window
    ui.geometry( "770x220" )

    # destroy SSH button
    serverButton.destroy()
    
    # state variables
    state_vars = []
    for mykey in stBounds:
        state_vars.append( var( mykey, real=True ) )

    # input variables
    inp_vars = []
    for mykey in inpBounds:
        inp_vars.append( var( mykey, real=True ) )

    # parameters
    params = []
    for mykey in parEqs:
        params.append( var( mykey, real=True ) )

    # parametric equation list
    eqlist = []
    for mykey in parEqs:
        eqlist.append( eval( mykey + "- nsimplify(" + parEqs[ mykey ] + ")" ) )

    # values of parameters and assignments
    valTuple = solve( eqlist, tuple(params) )
    if type( valTuple ) is list:
        vals =  list( valTuple[0] )
    elif type(valTuple) is dict:
        vals = list( valTuple.values() )
    asgn = []
    ind = 0
    for mypar in params:
        myval1, myval2 = fraction( vals[ ind ] )
        asgn.append( str( mypar ) + " = " + "Interval(" + str( myval1 ) + "," + str( myval1 ) +
                     ")/Interval(" + str( myval2 ) + "," + str( myval2 ) + ")" + ";" )
        ind += 1

    # dynamics in symbolic form
    dynEqs = {}
    for mykey in dynamics:
        dynEqs[ mykey ] =  eval( "nsimplify(" + dynamics[ mykey ] + ")" )

    # write state bounds
    lbobj = open( "src/pywrite/stlb.txt", "w" )
    ubobj = open( "src/pywrite/stub.txt", "w" )
    for mykey in stBounds:
        lbobj.write( str( stBounds[mykey][0] ) + " ")
        ubobj.write( str( stBounds[mykey][1] ) + " ")
    lbobj.close()
    ubobj.close()

    # write input bounds
    lbobj = open( "src/pywrite/inplb.txt", "w" )
    ubobj = open( "src/pywrite/inpub.txt", "w" )
    for mykey in inpBounds:
        lbobj.write( str( inpBounds[mykey][0] ) + " ")
        ubobj.write( str( inpBounds[mykey][1] ) + " ")
    lbobj.close()
    ubobj.close()    

    # process and write to c++ format
    writetocpp(dynEqs, state_vars, inp_vars, params, asgn, vals)

    # change edit state button
    editStButton[ "command" ] = editstate
    editInpButton[ "command" ] = editinp
    editDynamicsButton[ "command" ] = editdyn
    editParButton[ "command" ] = editpar

    processButton[ "text" ] = "Simulate"
    processButton[ "command" ] = simulate

    # place entry fields for simulation hyperparameters
    maxTimeEntry.place( x = startPos[0] + 500, y = startPos[1] )
    maxTimeEntryLabel.place( x = startPos[0] + 580, y = startPos[1] )
    logDivsEntry.place( x = startPos[0] + 500, y = startPos[1] + yGap )
    logDivsEntryLabel.place( x = startPos[0] + 580, y = startPos[1] + yGap )
    timeStepEntry.place( x = startPos[0] + 500, y = startPos[1] + 2*yGap )
    timeStepEntryLabel.place( x = startPos[0] + 580, y = startPos[1] + 2*yGap )

    # destroy load button
    loadButton.destroy()

    # compile c++ file
    if server[ "state" ] == "no":
        os.system( "make compile10" )
    elif server[ "state" ] == "yes":
        pth = server[ "path" ]
        # remove pywrite and results directory in host
        execstr = "ssh " + server[ "suffix" ] + " \"" + "rm -rf "
        execstr += pth + "NonlinearReachability/src/pywrite "
        execstr += pth + "NonlinearReachability/results; "

        # make fresh pywrite directory
        execstr += "mkdir "
        execstr += pth + "NonlinearReachability/src/pywrite; "

        # make fresh results directory
        execstr += "mkdir "
        execstr += pth + "NonlinearReachability/results\" "
        os.system( execstr )
        
        # copy files to pywrite directory
        execstr = "scp " + "src/pywrite/* " + server[ "suffix" ]
        execstr += ":" + pth + "NonlinearReachability/src/pywrite/"
        os.system( execstr )

        # compile inside server
        execstr = "ssh " + server[ "suffix" ] + " \"cd " + pth 
        execstr += "NonlinearReachability/;" + " make compile9; cd\""
        os.system( execstr )

processButton = tk.Button( ui, text = "Process", command = process,
                          font = tkFont.Font( size = 23 ),
                          highlightbackground = "green",
                          bg = "white" )
processButton.place( x = startPos[0] + 200, y = startPos[1]  )

#----------------------------------------------------------------------
# main tkinter loop
#----------------------------------------------------------------------
ui.mainloop()

