import sys
from PyQt5.QtWidgets import *
 
 
class Window(QWidget):
    right_uname = "Admin"
    right_pword = "Password1"
 
    def __init__(self):
        super(Window,self).__init__()
 
        self.init_ui()
 
    def init_ui(self):
        self.lbl_intro = QLabel('Welcome, please login')
        self.lbl_enter_username = QLabel('Username:')
        self.lbl_enter_password = QLabel('Password:')
        self.txt_enter_username = QLineEdit()
        self.txt_enter_password = QLineEdit()
        self.cb_login = QCheckBox('Stay logged in?')
        self.btn_login = QPushButton('Login')
 
 
        self.grid = QGridLayout()
        self.grid.setSpacing(5)
 
        self.grid.addWidget(self.lbl_intro, 1, 1)
 
        self.grid.addWidget(self.lbl_enter_username, 2, 0)
        self.grid.addWidget(self.txt_enter_username, 2, 1)
 
        self.grid.addWidget(self.lbl_enter_password, 3, 0)
        self.grid.addWidget(self.txt_enter_password, 3, 1)
 
        self.grid.addWidget(self.cb_login, 4, 1)
        self.grid.addWidget(self.btn_login, 5, 1)
 
 
        self.v_box = QVBoxLayout()
        self.v_box.addStretch(0)
        self.v_box.addLayout(self.grid)
        self.v_box.addStretch(0)
 
        self.h_box = QHBoxLayout()
        self.h_box.addStretch(0)
        self.h_box.addLayout(self.v_box)
        self.h_box.addStretch(0)
 
        self.setLayout(self.h_box)
 
        self.btn_login.clicked.connect(lambda: self.btn_login_clk(self.txt_enter_username, self.txt_enter_password, self.cb_login.isChecked(), self.lbl_intro))
 
 
        self.setWindowTitle('Login test')
 
        self.show()
 
    def btn_login_clk(self, username, password, cb, intro):
        if username.text() == self.right_uname and password.text() == self.right_pword:
            if cb:
                intro.setText('Welcome,' + ' ' + self.right_uname + ' ' + 'cb ticked')
            else:
                self.mw = MainWindow()
                self.hide()
                self.mw.show()
        else:
            intro.setText('Wrong username or password')
            self.clear_box()
 
    def clear_box(self):
        self.txt_enter_username.clear()
        self.txt_enter_password.clear()
        self.txt_enter_username.setFocus()
 
 
class MainWindow(Window):
 
    def __init__(self):
        super(MainWindow, self).__init__()
 
        self.init_ui()
 
    def init_ui(self):
        self.lbl_intro = QLabel('Main Window')
        self.lbl_user_logged = QLabel('Welcome,' + ' ' + self.right_uname)
        self.lbl_append = QLabel('Write something')
        self.txt_write_box = QLineEdit()
        self.btn_append = QPushButton('Append')
        self.btn_logout = QPushButton('Logout')
 
        layout = QVBoxLayout()
        layout.addWidget(self.lbl_intro)
        layout.addWidget(self.lbl_user_logged)
        layout.addWidget(self.lbl_append)
        layout.addWidget(self.txt_write_box)
        layout.addWidget(self.btn_append)
        layout.addWidget(self.btn_logout)
 
        self.setLayout(layout)
        self.setWindowTitle('Main')
 
        self.btn_append.clicked.connect(self.append_clk)
        self.btn_logout.clicked.connect(self.logout_action)
 
        self.show()
 
    def append_clk(self):
        self.appendee = open('text.txt', 'a')
        self.appendee.write('hry')
        self.appendee.close()
 
    def logout_action(self):
        self.close()
        a_window.show()
        a_window.clear_box()
 
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    a_window = Window()
    sys.exit(app.exec_())