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
        date +"%Y-%m-%d" > TODAY
        echo "v0.5.0" > WORKFLOW_VERSION
    >>>

    output {
        String analysis_date = read_string("TODAY")
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

task task_version_capture {
    input {
        Array[VersionInfo] version_array
        String workflow_name
        String workflow_version_path
        String project_name
        String analysis_date
        File version_capture_py
    }

    VersionInfoArray versions = { "versions": version_array }

    command <<<
        python ~{version_capture_py} \
            --versions_json ~{write_json(versions)} \
            --workflow_name ~{workflow_name} \
            --workflow_version ~{workflow_version_path} \
            --project_name ~{project_name} \
            --analysis_date ~{analysis_date}
    >>>

    output {
        File version_capture_file = "version_capture_${workflow_name}_${project_name}_${workflow_version_path}.csv"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "mchether/py3-bio:v4"
    }
}