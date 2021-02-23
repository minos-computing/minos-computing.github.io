#!/bin/bash
export SC_SIGNAL_WRITE_CHECK=DISABLE
NVDLA_EMU_BASE=/usr/local/nvdla
TUTORIAL_DIR=$PWD

num_instances=$1
if [ -z "$num_instances" ]; then
    num_instances=1
fi

mkdir -p nvdla_emulator_configs
base_port=6000

for i in `seq 0 $(( num_instances-1 ))`; do
port=$(( base_port + i ))
cat <<EOT > nvdla_emulator_configs/port${port}.lua
CPU = {
    library = "libqbox-nvdla.so",
    extra_arguments = "-machine virt -cpu cortex-a57 -machine type=virt -nographic -smp 1 -m 1024 -kernel Image --append \"root=/dev/vda\" -drive file=rootfs.ext4,if=none,format=raw,id=hd0 -device virtio-blk-device,drive=hd0 -fsdev local,id=r,path=.,security_model=none -device virtio-9p-device,fsdev=r,mount_tag=r -netdev user,id=user0,hostfwd=tcp::${port}-:6667, -device virtio-net-device,netdev=user0"
}

ram = {
    size = 1048576,
    target_port = {
        base_addr = 0xc0000000,
        high_addr = 0xffffffff
    }
}

nvdla = {
    irq_number = 176,
    csb_port = {
        base_addr = 0x10200000,
        high_addr = 0x1021ffff
    }
}
EOT
cd $NVDLA_EMU_BASE
expect ${TUTORIAL_DIR}/run_nvdla_emulator.exp ${TUTORIAL_DIR}/nvdla_emulator_configs/port${port}.lua &
cd ${TUTORIAL_DIR}
done

wait