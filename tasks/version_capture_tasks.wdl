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
        echo "v0.0.0" > WORKFLOW_VERSION
    >>>
    output {
        String workflow_version = read_string("WORKFLOW_VERSION")
    }
    runtime {
        memory: "512 MiB"
        cpu: 1
        docker: "ubuntu:focal"
    }
}
