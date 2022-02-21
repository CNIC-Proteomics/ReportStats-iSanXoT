Dim oShell

Dim SRC_HOME
SRC_HOME = Wscript.Arguments.Item(0)

Dim R_EXE, LAUNCHR, OPT
R_EXE = """" & SRC_HOME & "\R-Portable\App\R-Portable\bin\Rscript.exe" & """"
OPT = "--no-save --no-environ --no-init-file --no-restore --no-Rconsole"
LAUNCHR = """" & SRC_HOME & "\Launcher\launch.R" & """"

Dim CMD
CMD = R_EXE & " " & OPT & " " & LAUNCHR & " " & """" & SRC_HOME & """"

Set oShell = WScript.CreateObject ("WSCript.shell")
oShell.run "CMD /C " & """" & CMD & """", 0, True

Set oShell = Nothing