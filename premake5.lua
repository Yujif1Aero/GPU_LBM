function copydir(src_dir, dst_dir, filter, single_dst_dir)
    filter = filter or "**"
    src_dir = src_dir .. "/"
    print('copy "' .. src_dir .. filter .. '" to "' .. dst_dir .. '".')
    dst_dir = dst_dir .. "/"
    local dir = path.rebase(".", path.getabsolute("."), src_dir) -- root dir, relative from src_dir

    os.chdir(src_dir) -- change current directory to src_dir
    local matches = os.matchfiles(filter)
    os.chdir(dir) -- change current directory back to root

    local counter = 0
    for k, v in ipairs(matches) do
        local target = iif(single_dst_dir, path.getname(v), v)
        --make sure, that directory exists or os.copyfile() fails
        os.mkdir(path.getdirectory(dst_dir .. target))
        if os.copyfile(src_dir .. v, dst_dir .. target) then
            counter = counter + 1
        end
    end

    if counter == #matches then
        print(counter .. " files copied.")
        return true
    else
        print("Error: " .. counter .. "/" .. #matches .. " files copied.")
        return nil
    end
end

solution "Lbm"
	configurations { "Debug", "Release", "RelWithDebInfo" }
    filter {"platforms:x64", "configurations:Debug or configurations:DebugGpu"}
      targetsuffix "64D"
      defines {"DEBUG"}
      symbols "On"
    filter {"platforms:x64", "configurations:Release or configurations:RelWithDebInfo"}
      targetsuffix "64"
      defines {"NDEBUG"}
      optimize "On"
    filter {"platforms:x64", "configurations:RelWithDebInfo"}
      symbols "On"
	configuration {}
	
	language "C++"
	flags { "NoMinimalRebuild"  }
	flags { "MultiProcessorCompile" }
   characterset("MBCS")

	if os.istarget("windows") then
		defines{ "WIN32" }
		buildoptions { "/wd4305 /wd4244 /wd4661 /wd4996"  }
	end
	if os.istarget("linux") then
		defines{ "__LINUX__" }
		links {"X11", "Xi", "Xxf86vm", "Xrandr", "pthread" }
		buildoptions{"-msse3", "-mssse3" }
		buildoptions{"-std=c++17"}
	end
	
	platforms {"x64"}
	cppdialect "C++17"

	filter {"x64", "Debug"}
		targetsuffix "64D"
	filter {"x64", "Release"}
		targetsuffix "64"
		optimize "Speed"
	filter {}
--	tahoeRoot = "../../network/"
--	includedirs { tahoeRoot, tahoeRoot.."contrib/include/" }
	
	project "lbmApp"
		kind "ConsoleApp"
		location "./build"
		includedirs { "./" }
		if os.istarget("windows") then
			links{ "version" }
		end
		if os.istarget("linux") then
		end

		files { "./src/**.h", "./src/**.cpp" } 

		filter {"x64", "Debug"}
			targetdir "./dist/debug"
		filter {"x64", "Release"}
			targetdir "./dist/release"

