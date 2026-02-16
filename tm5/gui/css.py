#
#-- starting attempts to control layout/appearance of GUI elements
#
#
#-- https://panel.holoviz.org/how_to/styling/apply_css.html
#
setup_stylesheet = """
:host {
  font-size: 1.25em;
  background: #B6E0F0;
/*  border-radius: 5px; */
/*  border: 1px black solid; */
}

:host(.setup-base) {
  border-radius: 5px;
  border: 1px black solid;
}

:host(.setup-tracer) {
  background: #C4EEFE;
}
"""
precomp_stylesheet = """
:host(.precomp-left) {
  font-size: 1.25em;
  background: #C1D8E0;
  border-radius: 5px;
  border: 1px black solid;
}

:host(.precomp-right) {
  font-size: 1.25em;
  background: #D0EAF2;
  border-radius: 5px;
  border: 1px black solid;
}

:host(.green) .noUi-handle {
  background-color: green
}

:host(.blue) .noUi-handle {
  background-color: blue
}
"""

