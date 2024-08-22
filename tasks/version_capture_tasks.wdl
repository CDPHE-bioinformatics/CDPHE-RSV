version 1.0

struct VersionInfo {
    String software
    String docker
    String version
}

# workaround cromwell bug with read_json of Array
# https://github.com/openwdl/wdl/issues/409
struct VersionInfoArray {
    Array[VersionInfo] versions
}

task workflow_version_capture {
    input {
    }
    meta {
        description: "capture version release"
    }
    command <<<
        echo "v0.1.0" > WORKFLOW_VERSION
    >>>
    output {
        String workflow_version = read_string("WORKFLOW_VERSION")
        String workflow_version_path = sub(workflow_version, "\\.", "_")
    }
    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "ubuntu:focal"
    }
}
