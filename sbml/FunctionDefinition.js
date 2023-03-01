class FunctionDefinition {
    constructor( sbmlFunctionDefinition, libsbml) {
        this.sbmlFunctionDefinition = sbmlFunctionDefinition;
        this.functionName = this.sbmlFunctionDefinition.getName();
        this.id = this.sbmlFunctionDefinition.getId();
        this.argumentNames = []
        for (var i = 0; i < this.sbmlFunctionDefinition.getNumArguments(); i++) {
            this.argumentNames.push(this.sbmlFunctionDefinition.getArgument(i).getName());
        }
        this.body = libsbml.formulaToL3String(this.sbmlFunctionDefinition.getBody());
        this.libsbml = libsbml;
    }

    toString() {
        argumentCall = this.argumentNames.toString();
        callStr = `${this.id}(${argumentCall})`;
        return `${callStr}: ${this.body}`;
    }
}