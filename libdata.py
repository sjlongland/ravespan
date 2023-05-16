"""
#############################################################################
###########                    DATA VIEW                     ################
#############################################################################
"""
#################################################################
# This file is part of the RaveSpan program
# Copyright (C) 2017 Bogumil Pilecki
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################
#
# RaveSpan program includes the following files:
# * ravespan.py  - main RaveSpan GUI
# * librvc.py    - radial velocity curve plot
# * libspec.py   - spectrum related operations
# * libdata.py   - information about data points
# * libanal.py   - analysis window, velocity measurements
# * libdial.py   - various additional dialog windows
# * libcommpn.py - a bunch of small, useful functions 
# * rvspec       - executable python script to run RaveSpan
# * makelib      - custom installation bash script
# * work/specconv.py - a template program for spectra conversion to RaveSpan format
#
#################################################################



from PyQt5 import QtGui, QtCore, QtWidgets
from libcommon import *

class DATA_VIEW(QtWidgets.QMainWindow):
    """ Data View Window """
    def __init__(self, parent):
        super(DATA_VIEW, self).__init__(parent)
        self.setParent(parent)
        self.main = parent
        self.curve = self.main.curve
        self.spectrum = self.main.spectrum
        self.analysis = self.main.analysis

        self.setWindowTitle("Data View")
        self.main_frame = QtWidgets.QWidget() 
        
        self.Nv = self.main.Nv

        cmask0 = QtGui.QColor(233,233,233)
        cmask2 = QtGui.QColor(211,211,211)
        self.qcolors = {0: [cmask0,QtGui.QColor(222,222,255),cmask2],
                        1: [cmask0,QtGui.QColor(255,222,222),cmask2],
                        2: [cmask0,QtGui.QColor(222,255,222),cmask2],
                        3: [cmask0,QtGui.QColor(222,222,255),cmask2]}

        self.col_order = ["hjd", "rv1", "err1", "rv2", "err2", "rv3", "err3", "rv4", "err4", "instr", "snrat", "date", "w0", "wres","fname"]
        self.col_show={"hjd": True,
                      "rv1": True, "err1": True,
                      "rv2": True, "err2": True,
                      "rv3": False, "err3": False,
                      "rv4": False, "err4": False,
                      "instr": True,
                      "snrat": True,
                      "date": False,
                      "w0": False, "wres": False,
                      "fname": False}
                      
        self.yellow_cols = [self.col_order.index(x) for x in ["hjd", "instr", "snrat", "date", "w0", "wres", "fname"] if x in self.col_order]

        self.menubar = self.menuBar()
        self.mbMenu = {}
        self.mbMenu["main"] = self.menubar.addMenu('RV Curve')

        filemenu_order = ["open", "saveall", "saveselas", "--sep", "quit"]
        filemenu_items = {"open": ["Open", "Ctrl+O", 'Open object', lambda: 0],
            'saveall': ["Save all", "Ctrl+S", 'Save all comumns', self.curve.save_rvc_interactive],
            'saveselas': ["Save selected as...", None, 'Save only selected columns as...', self.save_selected_cols],
            'quit': ["Quit", "Ctrl+Q", 'Quit application', QtCore.SLOT('close()')]
        }

        for fmitem in filemenu_order:
            if fmitem == '--sep':
                self.mbMenu["main"].addSeparator()
            else:
                name, shortcut, statustip, eventfun = filemenu_items[fmitem]
                obj = QtWidgets.QAction(name, self)
                if shortcut is not None:
                    obj.setShortcut(shortcut)
                obj.setStatusTip(statustip)
                obj.triggered.connect(eventfun)
                self.mbMenu["main"].addAction(obj)

        self.mbMenu["partial"] = self.menubar.addMenu('Partial RVC')

        filemenu_order = ["residual", "orbit", "3rd_body", "puls1", "puls2"]
        filemenu_items = {"residual": ["Save residual", self.curve.save_residual],
                          "orbit": ["Save orbital", self.curve.save_orbital],
                          "3rd_body": ["Save third body", self.curve.save_3rd_body],
                          'puls1': ["Save pulsations A", self.curve.save_puls1],
                          'puls2': ["Save pulsations B", self.curve.save_puls2]
        }

        for fmitem in filemenu_order:
            if fmitem == '--sep':
                self.mbMenu["partial"].addSeparator()
            else:
                name, eventfun = filemenu_items[fmitem]
                obj = QtWidgets.QAction(name, self)
                obj.triggered.connect(eventfun)
                self.mbMenu["partial"].addAction(obj)

        self.morder = {}
        self.mitems = {}
        self.mobjects = {}
        self.mobjects["tools"]= {}
        self.mbMenu["tools"] = self.menubar.addMenu('Tools')

        self.morder["tools"] = ["err_rescale", 'calc_snrat', "show_vel1", "show_vel2", "show_vel3", "show_vel4",
                                "show_instr", "show_date", "show_wrange", "show_wres", "show_fname"]
        self.mitems["tools"] = {"err_rescale": ["Rescale errors", self.curve.rescale_errors, None],
                      "calc_snrat": ["Calculate S/N", self.curve.calculate_snrat, None],
                      "show_vel1": ["Show velocity 1", self.toggle_show_vel1, 1],
                      "show_vel2": ["Show velocity 2", self.toggle_show_vel2, 2],
                      "show_vel3": ["Show velocity 3", self.toggle_show_vel3, 0],
                      "show_vel4": ["Show velocity 4", self.toggle_show_vel4, 0],
                      "show_instr": ["Show instrument", self.toggle_show_instr, 1],
                      "show_date": ["Show date", self.toggle_show_date, 0],
                      "show_wrange": ["Show lamda_0", self.toggle_show_wrange, 0],
                      "show_wres": ["Show lamda res.", self.toggle_show_wres, 0],
                      "show_fname": ["Show file name", self.toggle_show_fname, 0],
        }

        for mitem in self.morder["tools"]:
            if mitem == '--sep':
                self.mbMenu["tools"].addSeparator()
            else:
                name, eventfun, checked = self.mitems["tools"][mitem]
                self.mobjects["tools"][mitem] = QtWidgets.QAction(name, self)
                self.mobjects["tools"][mitem].triggered.connect(eventfun)
                if checked is not None:
                    self.mobjects["tools"][mitem].setCheckable(True)
                    self.col_show["date"] = checked
                    self.mobjects["tools"][mitem].setChecked(checked)
                self.mbMenu["tools"].addAction(self.mobjects["tools"][mitem])

        self.scaling_factor_txt = self.menubar.addMenu("esf=%.3g"%self.curve.error_scaling)
        self.scaling_factor_txt.setEnabled(False)

        vbox = QtWidgets.QVBoxLayout()
        
        self.coldata = {"hjd":  ["HJD", 118, "%.5f"],
                        "rv1":  ["RV1",  75, "%.3f"],
                        "err1": ["ERR1", 60, "%.3f"],
                        "rv2":  ["RV2",  75, "%.3f"],
                        "err2": ["ERR2", 60, "%.3f"],
                        "rv3":  ["RV3",  75, "%.3f"],
                        "err3": ["ERR3", 60, "%.3f"],
                        "rv4":  ["RV4",  75, "%.3f"],
                        "err4": ["ERR4", 60, "%.3f"],
                        "instr":["INSTR", 100, "%s"],
                        "snrat":["S/N",  60, "%.3g"],
                        "date": ["DATE", 140, "%s"],
                        "w0":   ["WAVE_0", 90, "%.3f"],
                        "wres": ["W_RES", 65, "%.3g"],
                        "fname":["FNAME", 420, "%s"]
                        }
                        
        self.def_row = 16
        self.nrow = self.def_row
        self.def_rowheight = 21
        
        self.table = self.initialize_table()
                
        self.table.cellDoubleClicked.connect(self.point_selected)
        self.table.cellClicked.connect(self.point_selected_click)
        
        vbox.addWidget(self.table)
        
        self.main_frame.setLayout(vbox)          
        self.setCentralWidget(self.main_frame)   

        cwidths = self.get_visible_widths()
        self.setMinimumSize(int(sum(cwidths)*1.02) + 55, self.def_rowheight*(self.def_row+1)+40+20)


    def initialize_table(self):
        """ Initialize data table, set labels, widths, heights, set hidden columns. """
        table = QtWidgets.QTableWidget(self.def_row, len(self.col_order))
        table.setHorizontalHeaderLabels(self.get_labels())
        for i,key in enumerate(self.col_order):
            cwidth = self.coldata[key][1]
            table.setColumnWidth(i, cwidth)
            if not self.col_show[key]:
                table.hideColumn(i)   
        for i in range(self.def_row):
            table.setRowHeight(i, self.def_rowheight)
        
        table.setSelectionMode(QtGui.QAbstractItemView.NoSelection)
        
        table.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        
        return table


    def update_visible_cols(self):
        for i,key in enumerate(self.col_order):
            if self.col_show[key]:
                self.table.showColumn(i)
            else:
                self.table.hideColumn(i)   

    def toggle_show_vel1(self):
        state = not self.col_show["rv1"]
        self.col_show["rv1"] = state
        self.col_show["err1"] = state
        self.mobjects["tools"]["show_vel1"].setChecked(state)
        self.update_visible_cols()

    def toggle_show_vel2(self):
        state = not self.col_show["rv2"]
        self.col_show["rv2"] = state
        self.col_show["err2"] = state
        self.mobjects["tools"]["show_vel2"].setChecked(state)
        self.update_visible_cols()

    def toggle_show_vel3(self):
        state = not self.col_show["rv3"]
        self.col_show["rv3"] = state
        self.col_show["err3"] = state
        self.mobjects["tools"]["show_vel3"].setChecked(state)
        self.update_visible_cols()

    def toggle_show_vel4(self):
        state = not self.col_show["rv4"]
        self.col_show["rv4"] = state
        self.col_show["err4"] = state
        self.mobjects["tools"]["show_vel4"].setChecked(state)
        self.update_visible_cols()

    def toggle_show_instr(self):
        self.col_show["instr"] = not self.col_show["instr"]
        self.mobjects["tools"]["show_instr"].setChecked(self.col_show["instr"])
        self.update_visible_cols()

    def toggle_show_date(self):
        self.col_show["date"] = not self.col_show["date"]
        self.mobjects["tools"]["show_date"].setChecked(self.col_show["date"])
        self.update_visible_cols()

    def toggle_show_wrange(self):
        state = not self.col_show["w0"]
        self.col_show["w0"] = state
        self.mobjects["tools"]["show_wrange"].setChecked(state)
        self.update_visible_cols()

    def toggle_show_wres(self):
        self.col_show["wres"] = not self.col_show["wres"]
        self.mobjects["tools"]["show_wres"].setChecked(self.col_show["wres"])
        self.update_visible_cols()

    def toggle_show_fname(self):
        self.col_show["fname"] = not self.col_show["fname"]
        self.mobjects["tools"]["show_fname"].setChecked(self.col_show["fname"])
        self.update_visible_cols()
    

    def get_visible_widths(self):
        return [self.coldata[cid][1] for cid,cval in self.col_show.items() if cval]
        

    def get_labels(self):
        return [self.coldata[cid][0] for cid in self.col_order]

    
    def get_formats(self, col_id = None):
        if col_id is None:
            return [self.coldata[cid][2] for cid in self.col_id]
        else:
            return self.coldata[col_id][2]


    def clear_table(self):
        self.table.clearContents()

    def create_table(self):
        hjds = self.main.curve.get_data_for_id("hjd")
        self.nrow = max(self.def_row, len(hjds))
        self.table.setRowCount(self.nrow)
        
        for i in range(self.nrow):
            self.table.setRowHeight(i,self.def_rowheight)
            for j,cid in enumerate(self.col_order):
                twi = QtWidgets.QTableWidgetItem("")
                self.table.setItem(i,j,twi)    

    def fill_table(self):
        self.clear_table()
        self.create_table()
        self.update_table()


    def update_rowvelo(self, rownum, k):
        colkeys = ["rv"+str(k+1), "err"+str(k+1)]
        data = self.main.curve.get_data_for_ids(colkeys, ipoint = rownum)
        if data is None:
            return None
        vmask = self.main.curve.get_data_for_id("mask"+str(k+1), ipoint = rownum)
        item_col = self.qcolors[k][vmask]
        
        for key,ditem in zip(colkeys,data):
            if ditem is not None:
                cid, colnum = self.get_coldata(key)
                fmt = self.get_formats(cid)
    
                cell = self.table.item(rownum, colnum)
                if vmask == 0:
                    cell.setText("")
                else:
                    cell.setText(fmt % ditem)
                cell.setBackgroundColor(item_col)

    def update_rowvelos(self, rownum):
        for i in range(2):
            self.update_rowvelo(rownum, i)

    update_point = update_rowvelo
    update_row = update_rowvelos
        
    def get_coldata(self, col):
        if type(col) is str:
            if col not in self.col_order:
                print("Warning.")
                return
            cid = col
            colnum = self.col_order.index(col)
        else:
            if col > len(self.col_order):
                print("Warning.")
                return
            cid = self.col_order[col]
            colnum = col
        return cid, colnum        

    def update_col(self, col):
        """ col - column number or name """
        cid, colnum = self.get_coldata(col)

        data = self.main.curve.get_data_for_id(cid)
        if data is None:
            return
        fmt = self.get_formats(cid)
        
        if cid[:2] in ["rv", "er"]:
            vmask = self.main.curve.get_data_for_id("mask"+cid[-1:])
            k = int(cid[-1:])-1
            masked = True
        else:
            masked = False

        if DEVELOP: print(func_name(), cid, colnum, fmt, data)

        for rownum, ditem in enumerate(data):
            cell = self.table.item(rownum, colnum)

            if masked and vmask[rownum] == 0:
                cell.setText("")
            else:
                cell.setText(fmt % ditem)

            if masked:
                item_col = self.qcolors[k][vmask[rownum]]
                cell.setBackgroundColor(item_col)
 

    def update_cols(self, *cols):
        if len(cols)==0:
            return
        if type(cols[0]) in [list, tuple]:
            cols = cols[0]
        if DEVELOP: print(func_name(), cols)
        for col in cols:
            self.update_col(col)

    def update_table(self):
        self.update_cols(self.col_order)

    def select_point(self, row):
        for i in range(self.nrow):
            if i == row:
                col = QtGui.QColor(255,255,0)
            else:
                col = QtGui.QColor(255,255,255)
            for j in self.yellow_cols:
                self.table.item(i,j).setBackgroundColor(col)

    def point_selected(self, row, col):
        if col in self.yellow_cols:
            self.curve.select_datapoint(row)
        elif col>0 and col<self.main.Nv*2+1:
            self.curve.toggle_datapoint(row, (col-1)/2, plot_data=True)

    def point_selected_click(self, row, col):
        if col == 0 or col > self.main.Nv*2:
            return
        mods = QtWidgets.QApplication.keyboardModifiers()
        if mods == QtCore.Qt.AltModifier:
            self.curve.toggle_datapoint(row, (col-1)/2, plot_data=True)
        elif mods == QtCore.Qt.ControlModifier:
            self.curve.delete_datapoint(row, (col-1)/2, plot_data=True)
        elif mods == QtCore.Qt.ShiftModifier and col < 5:
            self.curve.swap_velos(row, plot_data=True)


    def save_selected_cols(self):
        """ Save selected columns. """
        cols_to_save = []
        for idx in self.table.selectionModel().selectedColumns():
            col = idx.column()
            cid = self.col_order[col]
            if self.col_show[cid]:
                cols_to_save.append(cid)
        self.curve.save_rvc_interactive(cols_to_save)
