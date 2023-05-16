"""
#######################################
#   A SET OF DIFFERENT DIALOG BOXES   #
#######################################
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



from PyQt5 import QtGui, QtCore
from libcommon import *

"""
####################################
#   INSTRUMENT FILTER DIALOG BOX   #
####################################
"""

class InstrumentDialog(QtGui.QDialog):
    """ Instrument Dialog Window - manage visibility of data according to instrument used """
    def __init__(self, parent=None):
        super(InstrumentDialog, self).__init__(parent)
        self.main = parent
        
        self.setWindowTitle("Instruments")
        
        vlayout = QtGui.QVBoxLayout()
        
        self.def_row = 6
        clabels = ["instrument", 'A', "res.", "lambda_0", 'N', 'shift']
        cwidths = {"instrument": 120, 'A': 25, "res.": 60, "lambda_0": 75, 'N': 40, 'shift': 55}
        cw_sum = sum(cwidths.values())
        self.cind = {k:i for i,k in enumerate(clabels)}
        def_col = len(clabels)

        self.checkboxes = []
        self.instr_list = []
        
        self.instr_table = QtGui.QTableWidget(self.def_row, def_col, self)
        self.instr_table.setHorizontalHeaderLabels(clabels)
        for k,i in self.cind.items():
            self.instr_table.setColumnWidth(i, cwidths[k])
        
        vlayout.addWidget(self.instr_table)

        button_box = QtGui.QDialogButtonBox()
        close_button = button_box.addButton(QtGui.QDialogButtonBox.Close)
        self.connect(close_button, QtCore.SIGNAL('clicked()'), QtCore.SLOT('close()'))
        vlayout.addWidget(button_box)

        self.setLayout(vlayout)
        self.setMinimumSize(cw_sum+50, 30*(self.def_row+1)+30+30)
        

    def update_table(self, instr_list, d_res, d_f0, d_n):

        self.instr_table.clearContents()        
        for i in range(len(self.checkboxes)):
            del self.checkboxes[0]

        self.instr_list = instr_list
        self.instr_table.setRowCount(max(self.def_row, len(instr_list)))
        
        for i, instr in enumerate(instr_list):
            chkbox = QtGui.QCheckBox()
            chkbox.setChecked(1)
            self.connect(chkbox, QtCore.SIGNAL('stateChanged(int)'), self.main.curve.plot_data)
            self.checkboxes.append(chkbox)
            self.instr_table.setCellWidget(i, 1, chkbox)
                        
            txt = QtGui.QTableWidgetItem(instr)
            txt.setFlags(QtCore.Qt.ItemFlags(QtCore.Qt.ItemIsEnabled))
            self.instr_table.setItem(i, 0, txt)

            for col, dic, fmt in [[2,d_res,"%.3f"], [3,d_f0,"%.2f"], [4,d_n, "%d"]]:
                tobj = QtGui.QTableWidgetItem(fmt%dic[instr])
                tobj.setFlags(QtCore.Qt.ItemFlags(QtCore.Qt.ItemIsEnabled))
                tobj.setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignCenter)
                self.instr_table.setItem(i, col, tobj)
            
            for k,v in {'shift':'0.0'}.items():
                obj = QtGui.QTableWidgetItem(v)
                obj.setTextAlignment(QtCore.Qt.AlignCenter)
                self.instr_table.setItem(i, self.cind[k], obj)

    def collect_colvalues(self, col):
        values = {}
        for i, instr in enumerate(self.instr_list):
            tval = self.instr_table.item(i, self.cind[col]).text()
            if '.' in tval:
                val = float(tval)
            else:
                val = int(tval)
            values[instr] = val
        return values

    def collect_active(self):
        active_instr = []
        for chkbox, instr in zip(self.checkboxes, self.instr_list):
            if chkbox.isChecked():
                active_instr.append(instr)
        return active_instr

    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()
    


"""
#######################################
#   EXTRA FITTING PARAMS DIALOG BOX   #
#######################################
"""


class FittingDialog(QtGui.QDialog):
    """ Fitting Dialog Window - additional parameters to fit: third body, pulsations """
    def __init__(self, parent=None):
        super(FittingDialog, self).__init__(parent)
        self.main = parent
        
        self.setWindowTitle("Fitting")
        
        vlayout = QtGui.QVBoxLayout()
        
        self.rvmodifiers_cb = {}
        rvmod_layout, rvmod_grid = {}, {}

        self.rvmod_order = ["3rd_body", "puls1", "puls2"]
        self.rvmod_parorder = {} 

        self.rvmodifiers_cb["3rd_body"] = QtGui.QCheckBox("third body")
        self.connect(self.rvmodifiers_cb["3rd_body"], QtCore.SIGNAL('toggled(bool)'), self.main.param_changed)
        
        rvmod_grid["3rd_body"] = QtGui.QGridLayout()
        
        tb_type_label = QtGui.QLabel("Type:")
        self.tb_type = QtGui.QComboBox()
        self.tb_type.addItems(["12+3","1+3", "2+3"])
        self.connect(self.tb_type, QtCore.SIGNAL('currentIndexChanged(QString)'), self.tb_mode_change)
        rvmod_grid["3rd_body"].addWidget(tb_type_label,0,1)
        rvmod_grid["3rd_body"].addWidget(self.tb_type,0,2,1,2)
        
        self.tb_cb = {}
        self.tb_spin = {}
        self.rvmod_parorder["3rd_body"] = ['per', 't0', 'ecc', 'aop', 'k12', 'k3']
        self.plims_d = {}
        self.tb_pdef_d = {}

        for par in [['per',  "P:",   10., 0.001,    39999.,  0.1, 6], \
                   [  't0',  "T0", 5000., -9999.,   49999.,  1.0, 5], \
                   [ 'ecc', "ecc",    0.,    0.,      0.97, 0.01, 4], \
                   [ 'aop', "aop",    0.,  -0.1,   dpi+0.1,  0.1, 4], \
                   [ 'k12', "K12",    0., -500.,      500.,  0.5, 3], \
                   [  'k3',  "K3",    0., -500.,      500.,  0.5, 3]]:
            key, text, val, minval, maxval, step, decis = par
            self.tb_pdef_d[key] = val
            self.plims_d[key] = [minval, maxval]
            self.tb_cb[key] = QtGui.QCheckBox(text, self)
            self.tb_spin[key] = set_QSpin(self, val, rmin=minval, rmax=maxval, step=step,decimals=decis)
            self.connect(self.tb_spin[key], QtCore.SIGNAL('valueChanged(double)'), self.tb_param_changed)

        self.tb_cb['k12'].setChecked(True)

        for i, par in enumerate(self.rvmod_parorder["3rd_body"]):
            rvmod_grid["3rd_body"].addWidget(self.tb_cb[par],i+1,1)  
            rvmod_grid["3rd_body"].addWidget(self.tb_spin[par],i+1,2,1,2)


        self.puls_fourharms = [{},{}]
        
        self.puls_cb = [{},{}]
        self.puls_spin = [{},{}]
        self.rvmod_parorder["puls1"] = ['per', 't0', 'a0']
        self.rvmod_parorder["puls2"] = self.rvmod_parorder["puls1"]
        self.pulsnharm_spin = [None,None]
        self.puls_clear_butt = [None,None]

        self.puls_pdef_d = {}

        fun_puls_event = [self.puls1_param_changed, self.puls2_param_changed]

        for i, ipulsmod in enumerate(["puls1", "puls2"]):
            self.rvmodifiers_cb[ipulsmod] = QtGui.QCheckBox("pulsations - "+chr(65+i))
            self.connect(self.rvmodifiers_cb[ipulsmod], QtCore.SIGNAL('toggled(bool)'), self.main.param_changed)
            rvmod_grid[ipulsmod] = QtGui.QGridLayout()

            for par in [['per', "P:",   10., 0.001,  39999., 0.01, 8], \
                       [  't0', "T0", 5000., -9999., 49999., 0.01, 5], \
                       [  'a0', "a0",    0., -9.99,    9.99,  0.1, 3]]:
                key, text, val, minval, maxval, step, decis = par
                self.puls_pdef_d[key] = val
                self.puls_cb[i][key] = QtGui.QCheckBox(text, self)
                self.puls_spin[i][key] = set_QSpin(self, val, rmin=minval, rmax=maxval, \
                                                   step=step, decimals=decis)
                self.connect(self.puls_spin[i][key], QtCore.SIGNAL('valueChanged(double)'), fun_puls_event[i])
            self.plims_d["a0"] = [-9.99, 9.99]

            self.puls_cb[i]['t0'].setChecked(True)
            pulsnharm_lab = QtGui.QLabel('nharm:', self)
            self.pulsnharm_spin[i] = set_QIntSpin(self,2,rmin=0,rmax=9)
            self.pulsnharm_spin[i].setMaximumWidth(90)
            self.connect(self.pulsnharm_spin[i], QtCore.SIGNAL('valueChanged(int)'), fun_puls_event[i])
            self.puls_clear_butt[i] = QtGui.QPushButton("Clear")
            self.puls_clear_butt[i].setMaximumWidth(70)
            self.puls_clear_butt[i].setAutoDefault(False)

            for j,ppar in enumerate(self.rvmod_parorder["puls1"]):
                rvmod_grid[ipulsmod].addWidget(self.puls_cb[i][ppar],j,1)
                rvmod_grid[ipulsmod].addWidget(self.puls_spin[i][ppar],j,2,1,2)

            rvmod_grid[ipulsmod].addWidget(pulsnharm_lab,3,1)
            rvmod_grid[ipulsmod].addWidget(self.pulsnharm_spin[i],3,2)
            rvmod_grid[ipulsmod].addWidget(self.puls_clear_butt[i],3,3)

        self.connect(self.puls_clear_butt[0], QtCore.SIGNAL('clicked()'), self.clear_puls_coefficients_1)
        self.connect(self.puls_clear_butt[1], QtCore.SIGNAL('clicked()'), self.clear_puls_coefficients_2)
        
        
        for rvmod in self.rvmod_order:
            rvmod_layout[rvmod] = QtGui.QHBoxLayout()
            rvmod_layout[rvmod].addSpacing(20)
            rvmod_layout[rvmod].addLayout(rvmod_grid[rvmod])
            vlayout.addWidget(self.rvmodifiers_cb[rvmod])
            vlayout.addLayout(rvmod_layout[rvmod])

        button_box = QtGui.QDialogButtonBox()
        close_button = button_box.addButton(QtGui.QDialogButtonBox.Close)
        self.connect(close_button, QtCore.SIGNAL('clicked()'), QtCore.SLOT('close()'))
        close_button.setAutoDefault(False)
        vlayout.addWidget(button_box)

        self.setLayout(vlayout)

    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()

    def tb_mode_change(self):
        mode = str(self.tb_type.currentText()).split('+')[0]
        self.tb_cb['k12'].setText("K"+mode)
        self.main.param_changed()

    def tb_param_changed(self):
        if self.rvmodifiers_cb["3rd_body"].isChecked():
            self.main.param_changed()
    
    def puls1_param_changed(self):
        if self.rvmodifiers_cb["puls1"].isChecked():
            self.main.param_changed()

    def puls2_param_changed(self):
        if self.rvmodifiers_cb["puls2"].isChecked():
            self.main.param_changed()        
    
    def get_active_modifiers(self):
        actmod = []
        actmodmask = []
        for imod in self.rvmod_order:
            if self.rvmodifiers_cb[imod].isChecked():
                actmod += [imod]
                actmodmask += [True]
            else:
                actmodmask += [False]
        return actmod

    def get_puls_nharm(self, i):
        return self.pulsnharm_spin[i].value()

    def clear_fields(self):
        for checkb in list(self.rvmodifiers_cb.values()):
            checkb.setChecked(False)
        
        self.tb_type.setCurrentIndex(0)
        
        for key, val in self.tb_pdef_d.items():
            self.tb_spin[key].setValue(val)
        
        for i in range(2):
            self.pulsnharm_spin[i].setValue(2)
        
        for key, val in self.puls_pdef_d.items():
            self.puls_spin[0][key].setValue(val)
            self.puls_spin[1][key].setValue(val)
        
        self.clear_puls_coefficients()
        

    def clear_puls_coefficients_1(self):
        self.puls_fourharms[0] = {}
        self.puls_clear_butt[0].setText("Clear")
        self.main.param_changed()

    def clear_puls_coefficients_2(self):
        self.puls_fourharms[1] = {}
        self.puls_clear_butt[1].setText("Clear")
        self.main.param_changed()

    def clear_puls_coefficients(self):
        self.puls_fourharms = [{},{}]
        for i in range(2):
            self.puls_clear_butt[i].setText("Clear")
        self.main.param_changed()

    def get_enabled_periods(self):
        val = []
        aff = [False]*self.main.Nv
        if self.rvmodifiers_cb["3rd_body"].isChecked():
            val.append(self.tb_spin["per"].value()/2.)
            tbtype_txt = self.tb_type.currentText()
            for i,si in enumerate(["1","2","3"]):
                if si in tbtype_txt:
                    aff[i] = True
            
        for i, ipulsmod in enumerate(["puls1", "puls2"]):
            if self.rvmodifiers_cb[ipulsmod].isChecked():
                val.append(self.puls_spin[i]["per"].value())
                aff[i] = True
        if len(val)==0:
            return [], []
        else:
            return val, aff

    def write_rv_mod_data(self, fobj):
        act_mod = self.get_active_modifiers()
        
        for rv_mod in self.rvmod_order:
            fobj.write("\n[%s]\n"%rv_mod)
            fobj.write("enabled = %s\n"%(rv_mod in act_mod))
            if rv_mod == "3rd_body":
                fobj.write("type   = %s\n"%self.tb_type.currentText())
            elif rv_mod in ["puls1", "puls2"]:
                fobj.write("nharm  = %s\n"%self.get_puls_nharm(int(rv_mod[-1:])-1))
            parvals = self.get_parvalues_for(rv_mod)
            for key, val in parvals.items():
                fobj.write("%-6s = %s\n"%(key,val))
        

    def get_parvalues_for(self, mod, withmask=False):
        parvalues = {}
        limits = {}
        parmask = {}
        if mod == "3rd_body":
            for key in self.rvmod_parorder["3rd_body"]:
                parvalues[key] = self.tb_spin[key].value()
                limits[key] = [self.tb_spin[key].minimum, self.tb_spin[key].maximum]
                parmask[key] = self.tb_cb[key].isChecked()
        if mod in ["puls1", "puls2"]:
            num = int(mod[4])-1
            for key in self.rvmod_parorder["puls1"]:
                parvalues[key] = self.puls_spin[num][key].value()
                limits[key] = [self.puls_spin[num][key].minimum, self.puls_spin[num][key].maximum]
                parmask[key] = self.puls_cb[num][key].isChecked()
            
            key, val = "a1", 0.0
            if "a1" in list(self.puls_fourharms[num].keys()):
                val = self.puls_fourharms[num]["a1"]
            parvalues[key], limits[key], parmask[key] = val, [None, None], True
            
            for i in range(1, self.pulsnharm_spin[num].value()):
                for c in ["a","b"]:
                    key, val = "%s%d"%(c,i+1), 0.0
                    if key in list(self.puls_fourharms[num].keys()):
                        val = self.puls_fourharms[num][key]
                    parvalues[key], limits[key], parmask[key] = val, [None, None], True
        
        if withmask:
            return parvalues, parmask
        else:                
            return parvalues
    
    
    def get_parlims_for_pars(self, parlist):
        parlims = []
        parlimited = []
        puls_four_lim = self.main.preferences.get_puls_four_lim()
        for p in parlist:
            if p in list(self.plims_d.keys()):
                parlims += [self.plims_d[p]]
                parlimited += [[1,1]]
            else: 
                parlims += [[-puls_four_lim, puls_four_lim]] 
                parlimited += [[1,1,]]
            
        return parlims, parlimited

    
    def set_parvalues_for(self, mod, pardict):
        if mod == "3rd_body":
            for key in list(pardict.keys()):
                self.tb_spin[key].setValue(pardict[key])

        nharm = 0
        if mod in ["puls1", "puls2"]:
            num = int(mod[4])-1
            for key in list(pardict.keys()):
                if key in self.rvmod_parorder["puls1"]:
                    self.puls_spin[num][key].setValue(pardict[key])
                else:
                    self.puls_fourharms[num][key] = pardict[key]
                    iharm = int(key[1:])
                    if iharm > nharm: nharm = iharm
            if nharm > 0:
                self.puls_clear_butt[num].setText("Clear (%d)"%nharm)
                    
    
    def get_parmask_for(self, mod):
        parmask = []
        return parmask
    
    def get_rvmods_and_opts(self):
        rvmods = []
        for key in self.rvmod_order:
            if self.rvmodifiers_cb[key].isChecked():
                rvmods += [key]
    
        rvopts = {}
        rvopts['tb_mode'] = str(self.tb_type.currentText())
        
        return rvmods, rvopts
    
    
    def set_fields(self, data):
        """ similar to 'set_parvalues_for' but uses 'enabled' key and more general 'data' dict """
        for key,val in data["3rd_body"].items():
            if key == "enabled":
                self.rvmodifiers_cb["3rd_body"].blockSignals(True)
                self.rvmodifiers_cb["3rd_body"].setChecked(val)
                self.rvmodifiers_cb["3rd_body"].blockSignals(False)
            elif key == "type":
                ci = self.tb_type.findText(val)
                if ci >= 0:
                    self.tb_type.setCurrentIndex(ci)
            else:
                self.tb_spin[key].blockSignals(True)
                self.tb_spin[key].setValue(val)
                self.tb_spin[key].blockSignals(False)
            
        for num, rv_mod in enumerate(["puls1", "puls2"]):
            nharm = 0
            for key,val in data[rv_mod].items():
                if key == "enabled":
                    self.rvmodifiers_cb[rv_mod].blockSignals(True)
                    self.rvmodifiers_cb[rv_mod].setChecked(val)
                    self.rvmodifiers_cb[rv_mod].blockSignals(False)
                elif key == "nharm":
                    self.pulsnharm_spin[num].blockSignals(True)
                    self.pulsnharm_spin[num].setValue(val)
                    self.pulsnharm_spin[num].blockSignals(False)
                else:
                    if key in self.rvmod_parorder["puls1"]:
                        self.puls_spin[num][key].blockSignals(True)
                        self.puls_spin[num][key].setValue(val)
                        self.puls_spin[num][key].blockSignals(False)
                    else:
                        self.puls_fourharms[num][key] = val
                        iharm = int(key[1:])
                        if iharm > nharm: nharm = iharm
                if nharm > 0:
                    self.puls_clear_butt[num].setText("Clear (%d)"%nharm)
    



"""
####################################
#  TEMPLATE SELECTION DIALOG BOX   #
####################################
"""

class TemplateDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(TemplateDialog, self).__init__(parent)
        self.main = parent

        self.tnum = 0
        self.pars = ['teff', 'logg', 'feh']
        self.data = [{'rot':0},{'rot':0}]
        self.combos = {}
        self.i_spins = {}
        self.i_steps = {"teff": 50, "logg": 0.25, "feh": 0.05}
        self.i_decs = {"teff": 0, "logg": 2, "feh": 2}
        self.int_exist = dict(list(zip(self.pars, [True,True,True])))   
        self.colors = [QtGui.QColor(0,0,0),QtGui.QColor(255,0,0),QtGui.QColor(0,210,0)]

        self.get_data_from_config()

        vlayout = QtGui.QVBoxLayout()
        layout = QtGui.QGridLayout()
        huplayout = QtGui.QHBoxLayout()
        vbotlayout = QtGui.QHBoxLayout()

        tmpl_label = QtGui.QLabel("Template:")

        self.tmpl=[]
        self.tmpl += [QtGui.QRadioButton("&1")]
        self.tmpl += [QtGui.QRadioButton("&2")]
        self.connect(self.tmpl[0], QtCore.SIGNAL('toggled(bool)'), self.switch_template)
        self.connect(self.tmpl[1], QtCore.SIGNAL('toggled(bool)'), self.switch_template)

        button_group = QtGui.QButtonGroup()        
        button_group.addButton(self.tmpl[0])
        button_group.addButton(self.tmpl[1])
        
        huplayout.addWidget(tmpl_label)
        huplayout.addWidget(self.tmpl[0])
        huplayout.addWidget(self.tmpl[1])

        tlabel = QtGui.QLabel("<html>T<sub>eff</sub>:</html>")
        tlabel.setToolTip("effective temperature")
        glabel = QtGui.QLabel("<html>log <i>g</i>:</html>")
        glabel.setToolTip("logarithm of gravitation")
        mlabel = QtGui.QLabel("[Fe/H]:")
        mlabel.setToolTip("metalicity")
        rlabel = QtGui.QLabel("<html><i>v</i>sin(<i>i</i>):</html>")
        rlabel.setToolTip("rotation corrected for inclination")

        for p in ['teff', 'logg', 'feh']:
            self.combos[p] = QtGui.QComboBox()

        self.rspin = set_QSpin(self,0,rmin=0,rmax=99,step=5, suffix=" km/s", decimals=0)
        self.rspin.setToolTip("<html>Setting <i>v</i>sin(<i>i</i>) differenet from 0 enables convolution of the selected template with a rotational profile. This process may cause some delay.</html>")

        self.pb_add_interpolated = QtGui.QPushButton("+")
        self.pb_add_interpolated.setToolTip("<html>Will add interpolated template to the template list. The more parameters differ from the grid values the more time it will take to make interpolated spectrum.</html>")
        self.pb_add_interpolated.setMaximumWidth(60)
        self.pb_add_interpolated.setEnabled(False)
        self.connect(self.pb_add_interpolated, QtCore.SIGNAL('clicked()'), self.add_interpolated )
        
        self.fill_combos()

        self.connect(self.combos['teff'], QtCore.SIGNAL('activated(QString)'), lambda v: self.combo_event('teff',v))
        self.connect(self.combos['logg'], QtCore.SIGNAL('activated(QString)'), lambda v: self.combo_event('logg',v))
        self.connect(self.combos['feh'], QtCore.SIGNAL('activated(QString)'), lambda v: self.combo_event('feh',v))
        self.connect(self.rspin, QtCore.SIGNAL('valueChanged(double)'), self.rot_changed)

        self.iname = QtGui.QLineEdit()
        self.iname.setToolTip("Interpolated spectrum name.")
        self.connect(self.iname, QtCore.SIGNAL('textChanged(QString)'), self.correct_iname)

        button_box = QtGui.QDialogButtonBox()
        close_button = button_box.addButton(QtGui.QDialogButtonBox.Close)
        self.connect(close_button, QtCore.SIGNAL('clicked()'), QtCore.SLOT('close()'))
        self.apply_button = button_box.addButton(QtGui.QDialogButtonBox.Apply)
        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.apply_params)
        vbotlayout.addWidget(button_box)

        layout.addWidget(tlabel, 1, 1)
        layout.addWidget(glabel, 2, 1)
        layout.addWidget(mlabel, 3, 1)
        layout.addWidget(rlabel, 4, 1)
        layout.addWidget(self.combos['teff'], 1, 2)
        layout.addWidget(self.combos['logg'], 2, 2)
        layout.addWidget(self.combos['feh'], 3, 2)
        layout.addWidget(self.rspin, 4, 2)
        layout.addWidget(self.pb_add_interpolated,4,3)
        layout.addWidget(self.iname,5,1,5,3)

        for p in self.pars:
            if len(self.grid[p])>0:
                self.i_spins[p] = set_QSpin(self, self.grid[p][(len(self.grid[p])-1)//2], \
                                        rmin=min(self.grid[p]), rmax=max(self.grid[p]), \
                                        step=self.i_steps[p], decimals=self.i_decs[p])
            else:
                self.i_spins[p] = set_QSpin(self, 0.)
                self.i_spins[p].setEnabled(False)

        self.connect(self.i_spins['teff'], QtCore.SIGNAL('valueChanged(double)'), lambda v: self.ispin_event('teff',v))
        self.connect(self.i_spins['logg'], QtCore.SIGNAL('valueChanged(double)'), lambda v: self.ispin_event('logg',v))
        self.connect(self.i_spins['feh'], QtCore.SIGNAL('valueChanged(double)'), lambda v: self.ispin_event('feh',v))


        layout.addWidget(self.i_spins['teff'], 1, 3)
        layout.addWidget(self.i_spins['logg'], 2, 3)
        layout.addWidget(self.i_spins['feh'], 3, 3)

        vlayout.addLayout(huplayout)
        vlayout.addLayout(layout)
        vlayout.addLayout(vbotlayout)

        self.setLayout(vlayout)
        self.setWindowTitle("Template X")


    def get_data_from_config(self):
        self.grid = {}
        if not os.path.exists(self.main.slibdir + '/slib.cfg'):
            print("INFO: Template spectra library (optional) not found. You can still use templates from ./templates directory.")
            #print "INFO: file '%s/slib.cfg' not found."%self.main.slibdir
            self.w0, self.res = "4000", "0.02"
            for p in self.pars:
                self.grid[p] = array([])
            self.main.disable_tpl_buttons()
            return
        with open(self.main.slibdir + '/slib.cfg') as slib:
            for line in slib:
                if line[13]=='=':
                    exec(line)
        self.w0 = str(wavelength_0)
        self.res = str(resolution)
        self.grid["teff"] = array(teff_grid)
        self.grid["logg"] = array(logg_grid)
        self.grid["feh"] = array(feh_grid)
        for p in self.pars:
            self.grid[p].sort()


    def fill_combos(self):
        """ FILL COMBO WIDGETS """
        for key in self.pars:
            if len(self.grid[key])==0: continue
            strvals = list(map(str, self.grid[key]))
            self.combos[key].addItems(strvals)
            for i in range(2):
                self.data[i][key] = strvals[(len(strvals)-1)//2]

    def switch_template(self, val):
        if self.tmpl[0].isChecked():
            self.tnum=0
        else:
            self.tnum=1
        for p in ['teff', 'logg', 'feh']:
            if p not in self.data[self.tnum]: continue
            d = str(self.data[self.tnum][p])
            i = self.combos[p].findText(d)
            self.combos[p].setCurrentIndex(i)
        self.rspin.setValue(self.data[self.tnum]['rot'])
        self.setWindowTitle("Template %d"%(self.tnum+1))
        self.set_interp_name()

    def combo_event(self, key, value):
        self.data[self.tnum][key] = str(value)
        self.i_spins[key].setValue(float(value))
        self.apply_button.setEnabled(True)

    def rot_changed(self, value):
        self.data[self.tnum]["rot"]=value
        self.apply_button.setEnabled(True)

    def ispin_event(self, key, value):
        if value not in self.grid[key]:
            self.int_exist[key] = False
        else:
            self.int_exist[key] = True
        self.pb_add_interpolated.setEnabled(not all(self.int_exist.values()))

        for key in self.pars:
            color = self.colors[int(not self.int_exist[key])]
            self.i_spins[key].setStyleSheet("QWidget { color: %s }" % color.name())            

    def select_template(self, tnum=0):
        self.setWindowTitle("Template %d"%(tnum+1))
        self.tnum=tnum
        self.tmpl[self.tnum].setChecked(1)

    def correct_iname(self,name):
        txt = self.iname.text()
        if '_' in txt:
            i = self.iname.cursorPosition()
            self.iname.setText(txt.replace('_','-'))
            self.iname.setCursorPosition(i)

    def set_interp_name(self,name=None):
        if name is None:
            name = self.main.obj_id_field.text()
        iid = name+"-"+str(self.tnum+1)
        self.iname.setText(iid)

    def apply_params(self):
        self.apply_button.setEnabled(False)
        idata = self.data[self.tnum]
        rotation = idata['rot']
        tid = "T%s g%s m%s r%s"%tuple([idata[key] for key in self.pars] + [str(rotation)])
        fname = "_".join([idata[key] for key in self.pars] + [self.w0] + [self.res])
        
        self.main.set_paratemplate(tid, fname, rotation, self.tnum)
        if self.main.preferences.apply_closes():
            self.hide()
        
    def add_interpolated(self):
        if all(self.int_exist.values()): return     
        self.pb_add_interpolated.setEnabled(False)
        for key in self.pars:
            if not self.int_exist[key]:
                self.i_spins[key].setStyleSheet("QWidget { color: %s }" % self.colors[2].name())

        igrid, ikeys = self.get_interpolation_grid()
        iparval_dict = dict([[p,[self.i_spins[p].value(),p in ikeys]] for p in self.pars])
        print(iparval_dict)
        itpl_files = self.get_template_grid(igrid)
        print(igrid)
        iid, iname = self.create_interpolated(iparval_dict, igrid, itpl_files, ikeys)
        self.main.add_template(iname)
        self.main.set_template(iid, iname, self.tnum)

    def get_fixed_value(self, key):
        diff = self.grid[key] - self.i_spins[key].value()
        i = argmin(abs(diff))
        return self.grid[key][i:i+1]
    
    def get_interp_values(self, key):
        diff = self.grid[key] - self.i_spins[key].value()
        i = argmin(abs(diff))
        if diff[i] == 0: raise ValueError
        if diff[i]<0:
            return self.grid[key][i:i+2]
        else:
            return self.grid[key][i-1:i+1]

    def get_interpolation_grid(self):
        interp_keys = []
        all_values = {}
        for key in self.pars:
            if self.int_exist[key]:
                all_values[key] = self.get_fixed_value(key)
            else:
                interp_keys += [key]
                all_values[key] = self.get_interp_values(key)

        return all_values, interp_keys

    def get_template_grid(self, all_values):
        dtfiles = []
        k1,k2,k3 = self.pars
        for v3 in all_values[k3]:
            for v2 in all_values[k2]:
                for v1 in all_values[k1]:
                    strlist = list(map(str, [v1,v2,v3])) + [self.w0] + [self.res]
                    tname = "_".join(strlist)
                    tpath = self.main.slibdir + "/" + tname
                    dtfiles += [fromfile(tpath,dtype=float32)]
        return dtfiles

    def linear_2specint(self, value, vals, specs):
        w = vals[1] - vals[0]
        w1 = (vals[1] - value)/w
        w2 = 1. - w1
        print("X:", w1, w2)
        ispec = specs[0]*w1 + specs[1]*w2
        return ispec

    def bilinear_4specint(self, values, vals1, vals2, specs):
        wx = vals1[1] - vals1[0]
        wy = vals2[1] - vals2[0]
        w1x = (vals1[1] - values[0])/wx
        w2x = 1. - w1x
        w1y = (vals2[1] - values[1])/wy
        w2y = 1. - w1y
        print("X:", w1x, w2x, " Y:", w1y, w2y)
        ispec = specs[0]*w1x*w1y + specs[1]*w2x*w1y + \
                specs[2]*w1x*w2y + specs[3]*w2x*w2y
        return ispec

    def trilinear_8specint(self, values, vals1, vals2, vals3, specs):
        wx = vals1[1] - vals1[0]
        wy = vals2[1] - vals2[0]
        wz = vals3[1] - vals3[0]
        w1x = (vals1[1] - values[0])/wx
        w2x = 1. - w1x
        w1y = (vals2[1] - values[1])/wy
        w2y = 1. - w1y
        w1z = (vals3[1] - values[2])/wz
        w2z = 1. - w1z
        print("X:", w1x, w2x, " Y:", w1y, w2y," Z:", w1z, w2z)
        ispec = specs[0]*w1x*w1y*w1z + specs[1]*w2x*w1y*w1z + \
                specs[2]*w1x*w2y*w1z + specs[3]*w2x*w2y*w1z + \
                specs[4]*w1x*w1y*w2z + specs[5]*w2x*w1y*w2z + \
                specs[6]*w1x*w2y*w2z + specs[7]*w2x*w2y*w2z
        return ispec
        
    def create_interpolated(self, iparval_dict, igrid, itpls, ikeys, limbdark=0.6, step_size=250):
        """ Creates file with interpolated template in the user template directory."""

        n_ikeys = len(ikeys)

        ivalues = [iparval_dict[p][0] for p in ikeys]
        print(ivalues, ikeys)
        
        if n_ikeys == 1:
            ik1 = ikeys[0]
            ispec = self.linear_2specint(ivalues[0], igrid[ik1], itpls)
        elif n_ikeys == 2:
            ik1 = ikeys[0]
            ik2 = ikeys[1]
            ispec = self.bilinear_4specint(ivalues, igrid[ik1], igrid[ik2], itpls)
        elif n_ikeys == 3:
            ik1 = ikeys[0]
            ik2 = ikeys[1]
            ik3 = ikeys[2]
            ispec = self.trilinear_8specint(ivalues, igrid[ik1], igrid[ik2], igrid[ik3], itpls)            

        vsini = self.data[self.tnum]["rot"]
        if vsini>0.0:
            res = float(self.res)
            w0  = float(self.w0)
            print(w0, res, vsini)
            ispec = convolve_with_rotational_profile(ispec, w0, res, vsini, limbdark, step_size)
        
        itid   = str(self.iname.text())
        itname = itid + "_%s_%s"%(self.w0, self.res)
        itpath =  self.main.tpl_dir + '/' + itname
        ispec.tofile(itpath)
        
        with open(itpath+".info", 'w') as finfo:
            for key in list(iparval_dict.keys()):
                finfo.write("%5s = %.2f\n"%(key,iparval_dict[key][0]))
        
        return itid, itname
            


"""
###################################
#     PREFERENCES DIALOG BOX      #
###################################
"""

class PreferencesDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(PreferencesDialog, self).__init__(parent)
        self.main = parent

        tabwidget = QtGui.QTabWidget(self)

        tab_order = ["basic", "adv"]
        tab_names = {"basic": "Basic", "adv": "Advanced"}

        tabs = {}
        tablayouts = {}
        for to in tab_order:
            tabs[to] = QtGui.QWidget()
            tablayouts[to] = QtGui.QVBoxLayout()
            tabs[to].setLayout(tablayouts[to])
            tabwidget.addTab(tabs[to], tab_names[to])
        
        self.group_order = ["gen", "tpl", "norm", "rvc", "anal", "spec", "fit",
                            "file", "advspec"]
        self.group_placement = {"adv": ["tpl", "norm", "advspec", "fit", "file", "--stretch"], 
                       "basic": ["gen", "rvc", "anal", "spec"] }
                       
        self.group_names = {"gen": "General:",
                       "tpl":  "Templates:",
                       "norm": "Normalization:",
                       "rvc":  "RV Curve (and fitting):",
                       "anal": "Analysis:",
                       "spec": "Spectrum:",
                       "fit": "Fitting:",
                       "file": "File options:",
                       "advspec": "Advanced spectrum:"
                       }
        self.groups = {}
        self.glayouts = {}
        self.gopts = {}
        for g in self.group_order:
            self.gopts[g] = {}

        for gid in self.group_order:
            self.groups[gid] = QtGui.QGroupBox(self.group_names[gid])
            self.glayouts[gid] = QtGui.QGridLayout(self.groups[gid])
        
        sect_items = {}
        
        sect_items["gen"] = {\
            "cb_update_plots_on_CA": [(1,1,1,1), 1, QtGui.QCheckBox, (1, "update plots on Calc All"),\
                                            ""] }
        
        sect_items["spec"] = {\
            "lab_min_res":     [(2,2,1,1), 1, QtGui.QLabel, ("<html>min. &lambda; step:</html>",),\
                                          "If data has higher resolution (lower step) it will be binned for calculations.\n0.0 means OFF."],
            "spin_min_res":    [(2,3,1,1), 1, set_QSpin, (0., 0.0, 1.0, 0.01), \
                                          "If data has higher resolution (lower step) it will be binned for calculations.\n0.0 means OFF."], 
            "lab_min_plotres": [(1,2,1,1), 1, QtGui.QLabel, ("<html>min. plot &Delta;&lambda;:</html>",),\
            "Minimum lambda step for spectrum plot.\nData is binned if it is lower (for faster plotting).\n0.0 means OFF."], \
            "spin_min_plotres":[(1,3,1,1), 1, set_QSpin, (0.1, 0.0, 1.0, 0.01), \
            "Minimum lambda step for spectrum plot.\nData is binned if it is lower (for faster plotting).\n0.0 means OFF."], \
            "cb_show_mask":    [(1,1,1,1), 1, QtGui.QCheckBox, (1, "show mask"),\
                                          "Show what part of spectrum will be used in the analysis."], \
            "cb_only_masked":  [(2,1,1,1), 1, QtGui.QCheckBox, (0, "only masked"),\
                                          ""],
            "cb_convert_onload":[(3,1,1,2), 0, QtGui.QCheckBox, (0, "convert res. on load"),\
                                          ""],
            "pb_update_spec":   [(3,3,1,1), 1, QtGui.QPushButton, ("Update",80), "Updates spectrum window."]}
        
        sect_items["rvc"] = {\
            "cb_autoswap":     [(1,1,1,1), 1, QtGui.QCheckBox, (1, "autoswap misidentified"),\
                                          "Autoswap data points for misidentified components."], \
            "cb_autodisable":  [(1,2,1,1), 1, QtGui.QCheckBox, (1, "autodisable outliers"),\
                                          "Automatically disable outlying points."], \
            "cb_show_err":     [(2,1,1,1), 1, QtGui.QCheckBox, (0, "show errors"),\
                                          "Show RV measurements errors."], \
            "cb_show_lc":      [(2,2,1,1), 1, QtGui.QCheckBox, (1, "show light curve"),\
                                          "Show light curve layer in RV Curve window."], \
            "lab_res_ran":     [(3,1,1,1), 1, QtGui.QLabel, ("res. view range",),\
                                          "Residuals view y-range (in km/s)."], \
            "spin_res_ran":    [(3,2,1,1), 1, set_QSpin, (0.81, 0.01, 9.9, 0.01), \
                                          "Residuals view y-range (in km/s)."] }

        sect_items["tpl"] = {\
            "lab_tpl_rot_limb":[(1,1,1,1), 1, QtGui.QLabel, ("limb darkening coefficient:",), \
                                        ""], \
            "spin_rot_limb":   [(1,2,1,1), 1, set_QSpin, (0.6, 0.0, 0.9, 0.1, 1), \
                                        ""], \
            "cb_apply_closes": [(2,1,1,2), 1, QtGui.QCheckBox, (1, "apply closes Template dialog box"),\
                                        "Autoswap data points for misidentified components."] }
        
        sect_items["norm"] = {
            "lab_binfill":[(1,1,1,1), 1, QtGui.QLabel, ("points per bin:",), ""],
            "ispin_binfill":   [(1,2,1,1), 1, set_QIntSpin, (1000, 100, 10000, 500), ""],
            "cb_norm_clean": [(1,3,1,1), 1, QtGui.QCheckBox, (1, "clean"), "Clean values higher than 1+noise for normalized spectrum."],
            "lab_mode_histbins":[(2,1,1,1), 1, QtGui.QLabel, ("bins for mode value:",), ""],
            "ispin_mode_histbins":   [(2,2,1,1), 1, set_QIntSpin, (25, 10, 100, 5), ""],
            "cb_norm_at": [(3,1,1,1), 1, QtGui.QCheckBox, (0, "normalize at:"), ""],
            "spin_norm_at":   [(3,2,1,1), 1, set_QSpin, (0.9, 0.25, 1.000, 0.005,3), ""],
            "pb_update_spec":   [(3,3,1,1), 1, QtGui.QPushButton, ("Update",100), "Updates spectrum window."]}

        sect_items["anal"] = {
            "lab_auto_rv":   [(1,1,1,1), 1, QtGui.QLabel, ("auto-detect RV:",),
                                     ""],
            "ispin_auto_rv": [(1,2,1,1), 1, set_QIntSpin, (2, 1, 4, 1),
                                        "How many velocities to detect automatically."],
            "lab_v_range":[(2,1,1,1), 1, QtGui.QLabel, ("v-range mul/add:",), 
                                        ""], 
            "spin_vran_mul":   [(2,2,1,1), 1, set_QSpin, (1.5, 1.0, 3.0, 0.1, 1), 
                                        ""], 
            "ispin_vran_add":   [(2,3,1,1), 1, set_QIntSpin, (30, 0, 99, 1), 
                                        ""], 
            "cb_model_id":   [(5,1,1,1), 1, QtGui.QCheckBox, (0, "ID from model"),
                                              "Use model for component identification."], 
            "cb_model_v":    [(5,2,1,1), 1, QtGui.QCheckBox, (0, "V from model"),
                                          "Use model to look for velocities."],
            "lab_ccf_fitfunc":   [(3,1,1,1), 1, QtGui.QLabel, ("CCF fit function:",),
                                     ""],
            "combo_ccf_fitfunc": [(3,2,1,2), 1, QtGui.QComboBox, (0, ["4th-order polynomial","gaussian"]),
                                          ""],
            "lab_bf_fitfunc":    [(4,1,1,1), 1, QtGui.QLabel, ("BF fit function:",),
                                     ""],
            "combo_bf_fitfunc": [(4,2,1,2), 1, QtGui.QComboBox, (1, ["rotational profile","gaussian","4th-order polynomial"]), ""]}

        sect_items["fit"] = {
            "lab_puls_four_lim":[(1,1,1,1), 1, QtGui.QLabel, ("Fourier coeff. limit (+/-):",),
                                        ""],
            "spin_fourier_lim":   [(1,2,1,1), 1, set_QSpin, (29.99, 0.01, 99.99, 0.1, 2),
                                        ""]}

        sect_items["file"] = {
           "lab_savespec_ftype": [(1,1,1,1), 1, QtGui.QLabel, ("Saved spectra file type:",), ""],
           "combo_savespec_ftype": [(1,2,1,2), 1, QtGui.QComboBox, (0, ["text", "binary"]), ""] }

        sect_items["advspec"] = {
            "lab_nsingular":[(1,1,1,1), 1, QtGui.QLabel, ("Max. rank of SVD:",), ""],
            "ispin_nsingular":   [(1,2,1,1), 1, set_QIntSpin, (501, 301, 1001, 100), "Maximum rank of the SVD matrices."],
            "lab_nkeepsvd":[(2,1,1,1), 1, QtGui.QLabel, ("# of SVD to keep:",), ""],
            "ispin_nkeepsvd":   [(2,2,1,1), 1, set_QIntSpin, (2, 1, 5, 1), "Number of SVD matrices to keep in memory."],
            }


        self.fitfuncs = {}
        self.fitfuncs["ccf"] = ["poly", "gauss"]
        self.fitfuncs["bf"] = ["rota", "gauss", "poly"]
        
                                        
        for section in list(sect_items.keys()):
            for item,values in sect_items[section].items():
                (x,y,w,h),enabled,widget,wpars,tooltip = values
                if item[:2] == "cb":
                    checked, name = wpars
                    self.gopts[section][item] = widget(name)
                    self.gopts[section][item].setChecked(checked)
                elif item[:2] in ["sp", "is"]:
                    if len(wpars) == 4:
                        dval, rmin, rmax, step = wpars
                        self.gopts[section][item] =  widget(self,dval,rmin=rmin,rmax=rmax,step=step)
                    elif len(wpars) == 5:
                        dval, rmin, rmax, step, deci = wpars
                        self.gopts[section][item] =  widget(self,dval,rmin=rmin,rmax=rmax,step=step,decimals=deci)
                elif item[:2] == "co":
                    cindex, citems = wpars
                    self.gopts[section][item] = widget()
                    self.gopts[section][item].addItems(citems)   
                    self.gopts[section][item].setCurrentIndex(cindex)
                elif item[:2] == "pb":
                    text, width = wpars
                    self.gopts[section][item] = widget(text)
                    self.gopts[section][item].setMaximumWidth(width)
                else:
                    name, = wpars
                    self.gopts[section][item] = widget(name)
                if tooltip:
                    self.gopts[section][item].setToolTip(tooltip)
                self.gopts[section][item].setEnabled(enabled)
                self.glayouts[section].addWidget(self.gopts[section][item],x,y,w,h)


        self.connect(self.gopts["gen"]["cb_update_plots_on_CA"], QtCore.SIGNAL('stateChanged(int)'), self.ca_update_changed)
        self.connect(self.gopts["spec"]["cb_show_mask"], QtCore.SIGNAL('stateChanged(int)'), self.set_mask_visibility)
        self.connect(self.gopts["spec"]["cb_only_masked"], QtCore.SIGNAL('stateChanged(int)'), self.main.spectrum.plot_data)
        self.connect(self.gopts["rvc"]["cb_show_err"], QtCore.SIGNAL('stateChanged(int)'), self.set_err_visibility)
        self.connect(self.gopts["anal"]["cb_model_id"], QtCore.SIGNAL('stateChanged(int)'), self.sync_model_v)
        self.connect(self.gopts["anal"]["cb_model_v"], QtCore.SIGNAL('stateChanged(int)'), self.sync_model_id)
        self.connect(self.gopts["spec"]["pb_update_spec"], QtCore.SIGNAL('clicked()'), self.main.spectrum.plot_data)
        self.connect(self.gopts["norm"]["pb_update_spec"], QtCore.SIGNAL('clicked()'), self.update_norm_spec)
        for field in ["ispin_binfill", "ispin_mode_histbins"]:
            self.connect(self.gopts["norm"][field], QtCore.SIGNAL('valueChanged(int)'), self.clear_db_norm)
        self.connect(self.gopts["norm"]["cb_norm_at"], QtCore.SIGNAL('stateChanged(int)'), self.clear_db_norm)
        self.connect(self.gopts["norm"]["spin_norm_at"], QtCore.SIGNAL('valueChanged(double)'), self.clear_db_norm_at)
        self.connect(self.gopts["advspec"]["ispin_nkeepsvd"], QtCore.SIGNAL('valueChanged(int)'), self.update_specdb_svd_max)

        for keytab,groups in self.group_placement.items():
            for gid in groups:
                if gid == "--stretch":
                    tablayouts[keytab].addStretch()
                else:
                    tablayouts[keytab].addWidget(self.groups[gid])

        mainlayout = QtGui.QVBoxLayout()
        mainlayout.addWidget(tabwidget)

        button_box = QtGui.QDialogButtonBox()
        save_button = button_box.addButton(QtGui.QDialogButtonBox.Save)
        close_button = button_box.addButton(QtGui.QDialogButtonBox.Close)
        self.connect(save_button, QtCore.SIGNAL('clicked()'), self.save_preferences)
        self.connect(close_button, QtCore.SIGNAL('clicked()'), QtCore.SLOT('close()'))
        mainlayout.addWidget(button_box)

        self.setLayout(mainlayout)
        self.setWindowTitle("Preferences")



    def update_main_window(self):
        self.ca_update_changed(self.gopts["gen"]["cb_update_plots_on_CA"].isChecked())
        self.main.curve.show_light_curve(self.gopts["rvc"]["cb_show_lc"].isChecked())

    def update_specdb_svd_max(self):
        """ Update SPECTRUM database maximum of SVD matrices to keep in memory. """
        self.main.spectrum.dbmax["svd"] = self.gopts["advspec"]["ispin_nkeepsvd"].value()
        #print self.main.spectrum.dbmax

    def sync_model_id(self, status):
        if status:
            self.gopts["anal"]["cb_model_id"].setChecked(True)
    
    def sync_model_v(self, status):
        if not status:
            self.gopts["anal"]["cb_model_v"].setChecked(False)

    def show_lc_changed(self, state):
        self.main.curve.show_light_curve(state, plot_data=True)

    def plot_rv_nostats(self):
        self.main.curve.plot_data(stats=False)

    def ca_update_changed(self,state):
        self.main.cb_ca_update.setChecked(state)
        
    def set_ca_update(self, state):
        self.gopts["gen"]["cb_update_plots_on_CA"].setChecked(state)

    def update_norm_spec(self):
        self.main.spectrum.plot_data()

    def clear_db_norm(self):
        self.main.spectrum.clear_db_norm()
    
    def clear_db_norm_at(self):
        if self.gopts["norm"]["cb_norm_at"].isChecked():
            self.main.spectrum.clear_db_norm()

    def get_res_range(self):
        return self.gopts["rvc"]["spin_res_ran"].value()

    def get_n_auto_vel(self):
        return self.gopts["anal"]["ispin_auto_rv"].value()

    def get_vrange_pars(self):
        mul = self.gopts["anal"]["spin_vran_mul"].value()
        add = self.gopts["anal"]["ispin_vran_add"].value()
        return mul, add

    def get_ID_from_model(self):
        return self.gopts["anal"]["cb_model_id"].isChecked()

    def get_v_from_model(self):
        return self.gopts["anal"]["cb_model_v"].isChecked()

    def get_useT2_for1D(self):
        return self.gopts["anal"]["cb_use_t2_in_1d"].isChecked()

    def get_minres(self):
        return self.gopts["spec"]["spin_min_res"].value()

    def get_plot_minres(self):
        return self.gopts["spec"]["spin_min_plotres"].value()


    def get_mask_only_view(self):
        return self.gopts["spec"]["cb_only_masked"].isChecked()

    def get_mask_visibility(self):
        return self.gopts["spec"]["cb_show_mask"].isChecked()

    def set_mask_visibility(self, value):
        self.main.spectrum.show_mask(self.gopts["spec"]["cb_show_mask"].isChecked())

    def set_err_visibility(self, value):
        self.main.curve.show_err(self.gopts["rvc"]["cb_show_err"].isChecked())

    def get_norm_params(self):
        binfill = self.gopts["norm"]["ispin_binfill"].value()
        modebins = self.gopts["norm"]["ispin_mode_histbins"].value()
        if self.gopts["norm"]["cb_norm_at"].isChecked():
            norm_at = self.gopts["norm"]["spin_norm_at"].value()
        else:
            norm_at = None 
        clean = self.gopts["norm"]["cb_norm_clean"].isChecked()
        return binfill, modebins, norm_at, clean

    def get_tpl_limbdark(self):
        return self.gopts["tpl"]["spin_rot_limb"].value()
    
    def apply_closes(self):
        return self.gopts["tpl"]["cb_apply_closes"].isChecked()

    def get_fit_function(self, method="ccf"):
        ind = self.gopts["anal"]["combo_%s_fitfunc"%method].currentIndex()
        return self.fitfuncs[method][ind]

    def get_num_singular_values(self):
        return self.gopts["advspec"]["ispin_nsingular"].value()

    def get_puls_four_lim(self):
        return self.gopts["fit"]["spin_fourier_lim"].value()

    def set_value(self, group, opt, val):
        if group in self.group_order and opt in self.gopts[group]:
            otype = opt.split('_')[0]
        else:
            return None
            
        try:
            if otype == "cb":
                self.gopts[group][opt].setChecked(eval(val))
            elif otype == "combo":
                self.gopts[group][opt].setCurrentIndex(int(val))
            elif otype == "spin":
                self.gopts[group][opt].setValue(float(val))
            elif otype == "ispin":
                self.gopts[group][opt].setValue(int(val))
        except:
            print("An error detected in the preferences file: ./preferences.")
            print("To fix it please go to the Preferences window and click the Save button.")

    def get_value(self, dopts, opt):
        otype = opt.split('_')[0]
        if otype == "cb":
            val = dopts[opt].isChecked()
        elif otype == "combo":
            val = dopts[opt].currentIndex()
        elif otype in ["spin", "ispin"]:
            val = dopts[opt].value()
        else:
            return None
        return str(val)

    def load_preferences(self):
        if not os.path.exists('preferences'):
            print("Warning: './preferences' file doesn't exist.")
            return
            
        with open("preferences","r") as fpref:
            for line in fpref:
                line = line.strip()
                if line == '':
                    continue
                if line[0]=='[' and line[-1]==']':
                    group = line[1:-1]
                elif '=' in line:
                    opt, eq, val = line.split()
                    self.set_value(group, opt, val)
        self.update_main_window()
        self.connect(self.gopts["rvc"]["cb_show_lc"], QtCore.SIGNAL('stateChanged(int)'), self.show_lc_changed)
        self.connect(self.gopts["rvc"]["spin_res_ran"], QtCore.SIGNAL('valueChanged(double)'), self.plot_rv_nostats)
        

    def save_preferences(self):
        with open("preferences","w") as fpref:
            for group in self.group_order:
                fpref.write("\n[%s]\n"%group)
                dopts = self.gopts[group]
                for opt in dopts:
                    value = self.get_value(dopts,opt)
                    if value is None:
                        continue
                    line = "%s = %s\n"%(opt,value)
                    fpref.write(line)
            fpref.write("\n")
        print("* preferences saved...")

    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()


"""
###################################
#       MASK EDIT DIALOG BOX      #
###################################
"""

class EditMaskDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(EditMaskDialog, self).__init__(parent)
        self.main = parent

        vlayout = QtGui.QVBoxLayout()

        self.editfield = QtGui.QPlainTextEdit()
        self.editfield.setMinimumWidth(400)
        self.editfield.setMinimumHeight(300)
        self.load_maskfile()
        vlayout.addWidget(self.editfield)

        button_box = QtGui.QDialogButtonBox()
        close_button = button_box.addButton(QtGui.QDialogButtonBox.Close)
        self.connect(close_button, QtCore.SIGNAL('clicked()'), QtCore.SLOT('close()'))
        self.apply_button = button_box.addButton(QtGui.QDialogButtonBox.Apply)
        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.apply_masks)
        vlayout.addWidget(button_box)

        self.setLayout(vlayout)
        self.setWindowTitle("Mask Edit")

    def load_maskfile(self):
        mask_fname = None
        for fname in ["masks.list", "masks.ranges", "spectrum.masks"]:
            if os.path.exists(self.main.datadir + "/" + fname):
                mask_fname = fname
                break
        if mask_fname is None:
            return
        with open(self.main.datadir + "/" + mask_fname) as smfile:
            text = smfile.read()
            self.editfield.setPlainText(text)

    def save_maskfile(self):
        with open(self.main.datadir + "/masks.list","w") as smfile:
            text = str(self.editfield.toPlainText())
            smfile.write(text.rstrip()+"\n")

    def apply_masks(self):
        self.save_maskfile()
        self.main.spectrum.load_masks()
        self.main.add_masks(clear=True)

    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()


"""
###################################
#       MASK EDIT DIALOG BOX      #
###################################
"""

class HelpDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(HelpDialog, self).__init__(parent)
        self.main = parent

        vlayout = QtGui.QVBoxLayout()

        self.txt = QtGui.QTextEdit()
        self.txt.setMinimumWidth(400)
        self.txt.setMinimumHeight(300)
        self.txt.setReadOnly(True)
        
        self.load_help()

        vlayout.addWidget(self.txt)

        button_box = QtGui.QDialogButtonBox()
        close_button = button_box.addButton(QtGui.QDialogButtonBox.Close)
        self.connect(close_button, QtCore.SIGNAL('clicked()'), QtCore.SLOT('close()'))
        vlayout.addWidget(button_box)

        self.setLayout(vlayout)
        self.setWindowTitle("Help")

    def load_help(self):
        if not os.path.exists("./HELP.html"):
            self.txt.setHtml("<h2 style='color: red'>Help file is missing !</h2>")
            return
        with open("HELP.html", "r") as hfile:
            helptxt = hfile.read()
            self.txt.setHtml(helptxt)

    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()
