# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(
    ['refitter.py'],
    datas=[
        ('../signaturesnet/data/data.xlsx', 'signaturesnet/data/'),
        ('../signaturesnet/data/mutation_type_order.xlsx', 'signaturesnet/data/'),
        ('../signaturesnet/data/datasets/example_input.csv', 'signaturesnet/data/datasets/'),
        
        ('../signaturesnet/trained_models/detector/init_args.json', 'signaturesnet/trained_models/detector/'),
        ('../signaturesnet/trained_models/detector/state_dict.zip', 'signaturesnet/trained_models/detector/'),
        ('../signaturesnet/trained_models/errorfinder/init_args.json', 'signaturesnet/trained_models/errorfinder/'),
        ('../signaturesnet/trained_models/errorfinder/state_dict', 'signaturesnet/trained_models/errorfinder/'),
        ('../signaturesnet/trained_models/finetuner_large/init_args.json', 'signaturesnet/trained_models/finetuner_large/'),
        ('../signaturesnet/trained_models/finetuner_large/state_dict.zip', 'signaturesnet/trained_models/finetuner_large/'),
        ('../signaturesnet/trained_models/finetuner_low/init_args.json', 'signaturesnet/trained_models/finetuner_low/'),
        ('../signaturesnet/trained_models/finetuner_low/state_dict.zip', 'signaturesnet/trained_models/finetuner_low/'),
        ('../signaturesnet/trained_models/generator/init_args.json', 'signaturesnet/trained_models/generator/'),
        ('../signaturesnet/trained_models/generator/state_dict', 'signaturesnet/trained_models/generator/'),
        
          ],
    pathex=[],
    binaries=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='refitter',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
