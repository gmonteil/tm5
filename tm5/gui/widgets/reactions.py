import panel as pn
import param
from typing import List


class ReactionSettings(pn.viewable.Viewer):
    """
    Block of widgets controlling *one* chemical reaction.
    Parameters:
        - reacname  ==> defined in the parent tracer class
        - shortname ==> defined in the parent tracer class
        - rate0 ==> user-defined? (not yet ...)
        - rate1 ==> user-defined? (not yet ...)
        - active ==> user defined (check-box)
        - domain ==> user defined (vertical domain)
        - field ==> user defined (from a list)
    """
    reacname = param.String(doc='name of the reaction (in the UI)')
    shortname = param.String(doc='name of the reaction (in the model)')
    rate0 = param.Number(doc="Reaction rate", label='rate a')
    rate1 = param.Number(doc='reaction rate', label='rate b')
    active = param.Boolean(default=True, doc='Enable or disable the reaction')
    domain = param.Selector(default='all', objects=['tropo', 'strato', 'all'], doc='Should the reaction be applied to the whole atmosphere or to a specific layer?')
    field = param.Selector(doc='version')

    def __init__(self, versions: List[str] = None, default_field: str = None, **param):
        super().__init__(**param)
        self.param.field.objects = versions
        self.param.field.doc = f'Version of the reactive component field to be used for the reaction {self.reacname}'
        self.field = versions[0]
        if default_field is not None: 
            self.field = default_field

    def __panel__(self):
        return pn.Row(
            pn.widgets.Checkbox.from_param(self.param.active, name=self.param.reacname),
            pn.widgets.Select.from_param(self.param.field) if len(self.param.field.objects) > 1 else None,
            #pn.widgets.FloatInput.from_param(self.param.rate0),
            #pn.widgets.FloatInput.from_param(self.param.rate1),
        )