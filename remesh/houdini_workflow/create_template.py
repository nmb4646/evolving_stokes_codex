"""Extract the approved green look from a legacy scene. Run with hython."""
import sys
from pathlib import Path
import hou
from PIL import Image

ROOT = Path(__file__).resolve().parent
source = Path(sys.argv[1]).resolve()
assets = ROOT / "assets"
assets.mkdir(exist_ok=True)
hou.hipFile.load(str(source), suppress_save_prompt=True, ignore_load_warnings=True)
hou.setUpdateMode(hou.updateMode.Manual)
obj = hou.node("/obj/geo1")
obj.setName("membrane")
obj.setWorldTransform(hou.Matrix4(1))
obj.node("measure1").setInput(0, obj.node("remesh1"))
for name in ("transform1", "transform2", "transform3"):
    node = obj.node(name)
    if node:
        node.destroy()
read = obj.node("python1")
read.setName("READ_MAT")
read.parm("python").set((ROOT / "sop_loader.py").read_text())
controls = obj.createNode("null", "CONTROLS")
group = hou.ParmTemplateGroup()
group.append(hou.StringParmTemplate("data_directory", "MAT Data Directory", 1, default_value=("",),
    string_type=hou.stringParmType.FileReference, file_type=hou.fileType.Directory))
group.append(hou.StringParmTemplate("source_indices", "Source Geometry Indices", 1, default_value=("0",)))
controls.setParmTemplateGroup(group)
controls.setColor(hou.Color((0.15, 0.65, 0.4)))
controls.setComment("Timeline frames 1..N read the listed geo indices from the data directory.")
controls.setGenericFlag(hou.nodeFlag.DisplayComment, True)
for name in ("Cut", "Wire"):
    node = hou.node("/mat/" + name)
    if node:
        node.destroy()
shader = hou.node("/mat/Membrane")
for parm in shader.parms():
    if parm.parmTemplate().type() == hou.parmTemplateType.String:
        try:
            if str(source.parent) in parm.evalAsString():
                parm.set("")
        except hou.Error:
            pass
Image.new("RGB", (8, 8), "white").save(assets / "white.png")
hou.node("/obj/cam1").parm("vm_background").set("$HIP/white.png")
rop = hou.node("/out/mantra_ipr")
rop.setName("render")
rop.parm("vm_picture").set("$HIP/frame_$F4.png")
rop.setParms({"vm_gamma": 1.0, "vm_usemaxthreads": 0, "vm_threadcount": 4, "soho_pipecmd": "mantra -j 4"})
obj.layoutChildren()
hou.setFrame(1)
hou.hipFile.save(str(assets / "green_membrane.hipnc"))
print("Saved", assets / "green_membrane.hipnc")
